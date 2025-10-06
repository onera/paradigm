/*============================================================================
 * Mesure des temps CPU et elapsed
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include "pdm_config.h"

#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <functional>
#include <iomanip>
#include <cstring>
#include <cmath>

#if defined (PDM_HAVE_GETRUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#elif defined(_POSIX_SOURCE)
#include <sys/times.h>
#include <unistd.h>
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_config.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure de mesure des temps d'execution
 *----------------------------------------------------------------------------*/

struct _pdm_timer_event_t {
  std::string event_name;
  std::string parent_name;
  std::string path_name;

  // Order list: Stores keys in the order they were first inserted (crucial for Call Tree sequence)
  std::vector<std::string> child_insertion_order;

  // Tree Structure: std::unordered_map for O(1) lookup and ownership
  std::unordered_map<std::string, std::unique_ptr<_pdm_timer_event_t>> children;

  // Node Statistics
  long   n_call          = 0;
  double t1              = 0.; // Start time of the active call
  double t_run_inclusive = 0.; // Total inclusive time
  double t_run_exclusive = 0.; // Exclusive time (calculated at dump)
  double t_children_sum  = 0.; // Sum of children's inclusive times
  double t_sync_entry    = 0.; // Synchronization wait time at the start (Entry Barrier Overhead)
  double t_sync_exit     = 0.; // Synchronization wait time at the end (Exit Barrier Overhead)

  // Recursion Management
  int    is_active_count = 0;

  _pdm_timer_event_t() = default;
  _pdm_timer_event_t(_pdm_timer_event_t const&) = delete;
  _pdm_timer_event_t& operator=(_pdm_timer_event_t const&) = delete;
};

  // double t_adim               = 0.;
struct _pdm_global_stat_t {
  long   n_call              = 0;
  double t_sum_inclusive     = 0.0;
  double t_sum_exclusive     = 0.0;
  double t_sum_sync_entry    = 0.0;
  double t_sum_sync_exit     = 0.0;

  double t_mean_run_inclusive = 0.0;
  double t_mean_run_exclusive = 0.0;
  double t_mean_sync_entry    = 0.0;
  double t_mean_sync_exit     = 0.0;

  // --- Min Inclusive
  double t_min_run_inclusive = HUGE_VAL;
  int    rank_min_inclusive  = -1;

  // --- Max Inclusive
  double t_max_run_inclusive = 0.0;
  int    rank_max_inclusive  = -1;

  // --- Min Exclusive
  double t_min_run_exclusive = HUGE_VAL;
  int    rank_min_exclusive  = -1;

  // --- Max Exclusive
  double t_max_run_exclusive = 0.0;
  int    rank_max_exclusive  = -1;

  // --- Min Sync Entry
  double t_min_sync_entry    = HUGE_VAL;
  int    rank_min_sync_entry = -1;

  // --- Max Sync Entry
  double t_max_sync_entry    = 0.0;
  int    rank_max_sync_entry = -1;

  // --- Min Sync Exit
  double t_min_sync_exit    = HUGE_VAL;
  int    rank_min_sync_exit = -1;

  // --- Max Sync Exit
  double t_max_sync_exit    = 0.0;
  int    rank_max_sync_exit = -1;
};


struct _pdm_timer_t {

  double  t_cpu;                   /* Temps CPU cumule */
  double  t_elapsed;               /* Temps elapsed cumule */
#if defined (PDM_HAVE_GETRUSAGE)
  double  t_cpu_u;                 /* Temps CPU utilisateur */
  double  t_cpu_s;                 /* Temps CPU system */
  double  t_cpu_debut;
  double  t_cpu_u_debut;
  double  t_cpu_s_debut;
#else
  clock_t t_cpu_debut;             /* Marque de debut de mesure
                                      du temps CPU */
#endif
  struct timeval t_elaps_debut;    /* Marque de debut de mesure
                                      du temps elapsed */
  int     indic;                   /* Indique si une mesure d'une tranche est en cours */


  PDM_MPI_Comm comm;
  _pdm_timer_event_t root_event;
  std::vector<_pdm_timer_event_t*> call_stack;
};

/*============================================================================
 * Definition des fonctions locales
 *============================================================================*/

/**
 * @brief Obtient le parent actif directement à partir de la pile.
 * @return Le pointeur vers l'événement parent, ou &timer->root_event si la pile est vide.
 */
static
_pdm_timer_event_t*
get_active_parent
(
  PDM_timer_t *timer
)
{
  if (timer->call_stack.empty()) {
    return &timer->root_event;
  }
  // The last element in the stack is the parent by construction
  return timer->call_stack.back();
}


/**
 * @brief Traverses the tree to calculate max width for the name column.
 */
static
void
calculate_max_widths(_pdm_timer_event_t& node, int current_indent_length, size_t& max_name_width) {
  if (node.event_name == "__ROOT__") {
    for (auto& child_name : node.child_insertion_order) {
      calculate_max_widths(*node.children.at(child_name), current_indent_length, max_name_width);
    }
    return;
  }

  // Indentation uses 3 chars per level ("  |")
  size_t name_length = node.event_name.length() + (size_t)current_indent_length;

  if (name_length > max_name_width) {
    max_name_width = name_length;
  }

  // Recurse
  for (auto& child_name : node.child_insertion_order) {
    calculate_max_widths(*node.children.at(child_name), current_indent_length + 3, max_name_width);
  }
}


/**
 * @brief Formats a single timer line (event statistics) into an aligned string.
 */
std::string format_timer_line(_pdm_timer_event_t& node,
                              const std::string& indented_name,
                              size_t name_width,
                              int time_width,
                              int ncall_width) {

    // 1. Calculate Exclusive Time
    double children_time_sum = 0.0;
    for (auto& child_name : node.child_insertion_order) {
      children_time_sum += node.children.at(child_name)->t_run_inclusive;
    }
    double t_exclusive = node.t_run_inclusive - children_time_sum;

    std::stringstream ss;

    // --- Column 1: Event Name (with Indentation) ---
    ss << std::left << std::setw(name_width) << indented_name;

    // --- Column 2: N_Call ---
    ss << std::right << std::setw(ncall_width) << node.n_call;

    // --- Column 3: T_Inclusive ---
    ss << " |"
       << std::right << std::setw(time_width)
       << std::fixed << std::setprecision(6)
       << node.t_run_inclusive;

    // --- Column 4: T_Exclusive ---
    ss << " |"
       << std::right << std::setw(time_width)
       << std::fixed << std::setprecision(6)
       << t_exclusive;

    // --- Column 5: T_Sync_Entry ---
    ss << " |"
       << std::right << std::setw(time_width)
       << std::fixed << std::setprecision(6)
       << node.t_sync_entry;

    // --- Column 6: T_Sync_Exit ---
    ss << " |"
       << std::right << std::setw(time_width)
       << std::fixed << std::setprecision(6)
       << node.t_sync_exit;

    return ss.str();
}

/**
 * @brief Traverses the tree recursively, formats lines, and collects them in a vector.
 */
void traverse_and_print_aligned(_pdm_timer_event_t& node,
                                int depth,
                                size_t name_width,
                                int time_width,
                                int ncall_width,
                                std::vector<std::string>& lines) {

  if (node.event_name == "__ROOT__") {
    for (auto& child_name : node.child_insertion_order) {
      traverse_and_print_aligned(*node.children.at(child_name), depth, name_width, time_width, ncall_width, lines);
    }
    return;
  }

  // 1. Format the name with indentation
  std::string indent = "";
  for (int i = 0; i < depth; ++i) {
    indent += "  |";
  }
  std::string indented_name = indent + node.event_name;

  // 2. Format the line using the dedicated function and store it
  lines.push_back(format_timer_line(node, indented_name, name_width, time_width, ncall_width));

  // 3. Recurse (using the insertion order)
  for (auto& child_name : node.child_insertion_order) {
    traverse_and_print_aligned(*node.children.at(child_name), depth + 1, name_width, time_width, ncall_width, lines);
  }
}

/**
 * @brief Collects all nodes into a single vector for flat mode printing.
 */
void collect_all_nodes(_pdm_timer_event_t* n, std::vector<_pdm_timer_event_t*>& all_nodes) {
  if (n->event_name != "__ROOT__") {
    all_nodes.push_back(n);
  }
  for (auto& child_name : n->child_insertion_order) {
    collect_all_nodes(n->children.at(child_name).get(), all_nodes);
  }
}


/**
 * @brief Recursively calculates exclusive time and dumps to JSON (using insertion order).
 */
void calculate_exclusive_and_dump(_pdm_timer_event_t& node, int indent, FILE* fp, bool& first_child) {

  if (node.event_name != "__ROOT__") {
    if (!first_child) {
      fprintf(fp, ",\n");
    }
    first_child = false;

    // --- Calculate Exclusive Time ---
    double children_time_sum = 0.0;
    for (auto& child_name : node.child_insertion_order) {
      children_time_sum += node.children.at(child_name)->t_run_inclusive;
    }
    node.t_run_exclusive = node.t_run_inclusive - children_time_sum;


    // --- JSON Serialization ---
    for (int i = 0; i < indent * 2; ++i) fprintf(fp, " ");
    fprintf(fp, "{\n");

    for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
    fprintf(fp, "\"name\": \"%s\",\n", node.event_name.c_str());

    for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
    fprintf(fp, "\"n_call\": %ld,\n", node.n_call);

    for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
    fprintf(fp, "\"t_inclusive\": %12.5e,\n", node.t_run_inclusive);

    for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
    fprintf(fp, "\"t_exclusive\": %12.5e,\n", node.t_run_exclusive);

    for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
    fprintf(fp, "\"t_sync_entry\": %12.5e,\n", node.t_sync_entry); // Entry

    for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
    fprintf(fp, "\"t_sync_exit\": %12.5e", node.t_sync_exit); // Exit

    // Add children if present
    if (!node.children.empty()) {
      fprintf(fp, ",\n");
      for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
      fprintf(fp, "\"children\": [\n");

      bool current_level_first_child = true;
      for (auto& child_name : node.child_insertion_order) {
        calculate_exclusive_and_dump(*node.children.at(child_name), indent + 2, fp, current_level_first_child);
      }
      fprintf(fp, "\n");
      for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
      fprintf(fp, "]");
    }

    fprintf(fp, "\n");
    for (int i = 0; i < indent * 2; ++i) fprintf(fp, " ");
    fprintf(fp, "}");
  }
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

PDM_timer_t*
PDM_timer_create2
(
  PDM_MPI_Comm comm
)
{
  PDM_timer_t *timer = new PDM_timer_t();

  timer->comm = comm;

  // Init root to begin stack
  timer->root_event.event_name  = "__ROOT__";
  timer->root_event.parent_name = "__NONE__";
  timer->root_event.path_name   = "";
  return timer;
}



void
PDM_timer_start
(
        PDM_timer_t *timer,
  const char        *name,
        int         force_synchro
)
{
  std::string current_name(name);
  _pdm_timer_event_t* parent_node = get_active_parent(timer);
  _pdm_timer_event_t* event_ptr;

  // 1. Get or create the node
  if (parent_node->children.count(current_name) == 0) {

    // Record the insertion order
    parent_node->child_insertion_order.push_back(current_name);

    // Allocation via unique_ptr
    parent_node->children[current_name] = std::make_unique<_pdm_timer_event_t>();
    event_ptr = parent_node->children.at(current_name).get();

    // Initialization
    event_ptr->event_name  = current_name;
    event_ptr->parent_name = parent_node->event_name;
    event_ptr->path_name   = parent_node->path_name + "/" + current_name;
  } else {
    // If the child exists, retrieve it.
    event_ptr = parent_node->children.at(current_name).get();
  }

  // 2. Measure and execute Entry Barrier
  if (force_synchro == 1) {
    double t_sync_start = PDM_MPI_Wtime();
    PDM_MPI_Barrier(timer->comm);
    // Accumulate the synchronization wait time at entry
    event_ptr->t_sync_entry += (PDM_MPI_Wtime() - t_sync_start);
  }

  // 3. Time and Stack Management
  _pdm_timer_event_t& event = *event_ptr;

  if (event.is_active_count == 0) {
    event.t1 = PDM_MPI_Wtime();
  }

  event.n_call++;
  event.is_active_count++;

  // Push the POINTER of the current event onto the stack
  timer->call_stack.push_back(event_ptr);
}


void
PDM_timer_end
(
        PDM_timer_t *timer,
  const char        *name,
        int          force_synchro
)
{
  std::string current_name(name);

  // 1. Récupérer le n?ud qui vient de se terminer (sommet de la pile)
  if (timer->call_stack.empty()) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_TIMER ERROR: PDM_timer_end(%s) called but the stack is empty \n", name);
    return;
  }

  _pdm_timer_event_t* current_node = timer->call_stack.back();
  if (current_node->event_name != current_name) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_TIMER ERROR: PDM_timer_end (%s) called but active timer is %s \n", name, current_node->event_name);
    return;
  }

  // Pop the finished event from the stack
  timer->call_stack.pop_back();

  // Retrieve the new parent node
  _pdm_timer_event_t* parent_node = get_active_parent(timer);

  // NEW: Measure and execute Exit Barrier
  if (force_synchro == 1) {
    double t_sync_start = PDM_MPI_Wtime();
    PDM_MPI_Barrier(timer->comm);
    // Accumulate the synchronization wait time at exit
    current_node->t_sync_exit += (PDM_MPI_Wtime() - t_sync_start);
  }

  current_node->is_active_count--;

  if (current_node->is_active_count == 0) {
    double dt = PDM_MPI_Wtime() - current_node->t1;
    current_node->t_run_inclusive += dt;

    // Update parent's children sum
    parent_node->t_children_sum += dt;
  }
}


/**
 * @brief Generates the full formatted report string (header + lines) for logging.
 * The returned string must be freed by the user using PDM_timer_free_string.
 * NOTE: The user is responsible for freeing the returned char* using the dedicated function.
 * * @param timer The timer instance.
 * @param mode 0: Hierarchical (Indented). 1: Flat/Raw (All nodes, non-indented).
 * @return Dynamically allocated C-string containing the full report.
 */
char* PDM_timer_get_report_string(PDM_timer_t *timer, int mode) {
  std::stringstream report_stream;
  if (timer->root_event.children.empty()) {
    report_stream << "PDM_TIMER: No events recorded.\n";
    // Return a copy of the stringstream content
    std::string temp_str = report_stream.str();
    char *cstr = new char[temp_str.length() + 1];
    std::strcpy(cstr, temp_str.c_str());
    return cstr;
  }

  // --- Phase 1: Calculate Dynamic Widths ---
  size_t max_name_width = 11;

  if (mode == 0) {
    calculate_max_widths(timer->root_event, 0, max_name_width);
  } else {
    std::vector<_pdm_timer_event_t*> all_nodes;
    collect_all_nodes(&timer->root_event, all_nodes);
    for(auto* node : all_nodes) {
      if (node->event_name.length() > max_name_width) {
        max_name_width = node->event_name.length();
      }
    }
  }

  const size_t MAX_COL_LIMIT = 80;
  max_name_width = std::min(max_name_width + 2, MAX_COL_LIMIT);

  const int N_CALL_COL_WIDTH = 10;
  const int TIME_COL_WIDTH = 18;
  const int TOTAL_WIDTH = max_name_width + N_CALL_COL_WIDTH + (TIME_COL_WIDTH * 4) + 8;

  // --- Phase 2: Build Header ---

  report_stream << "\n" << std::string(TOTAL_WIDTH, '=') << "\n";
  report_stream << "PDM TIMER REPORT (" << (mode == 1 ? "FLAT/RAW MODE" : "HIERARCHICAL MODE") << ")\n";
  report_stream << std::string(TOTAL_WIDTH, '=') << "\n";

  // Column Header
  report_stream << std::left << std::setw(max_name_width) << "Event Name (Path)";
  report_stream << std::right << std::setw(N_CALL_COL_WIDTH) << "Calls";
  report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "Time Inclusive (s)";
  report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "Time Exclusive (s)";
  report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "Entry Wait (s)";
  report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "Exit Wait (s)";
  report_stream << "\n" << std::string(TOTAL_WIDTH, '-') << "\n";

  // --- Phase 3: Build Content Lines ---

  std::vector<std::string> content_lines;

  if (mode == 0) { // Hierarchical Mode
    for (auto& child_name : timer->root_event.child_insertion_order) {
      traverse_and_print_aligned(*timer->root_event.children.at(child_name), 0, max_name_width, TIME_COL_WIDTH, N_CALL_COL_WIDTH, content_lines);
    }
  } else { // Flat/Raw Mode
    std::vector<_pdm_timer_event_t*> all_nodes;
    collect_all_nodes(&timer->root_event, all_nodes);

    for(auto* node : all_nodes) {
      // For flat mode, the indented name is just the event name (no prefix)
      content_lines.push_back(format_timer_line(*node, node->event_name, max_name_width, TIME_COL_WIDTH, N_CALL_COL_WIDTH));
    }
  }

  // Append all content lines
  for(const auto& line : content_lines) {
    report_stream << line << "\n";
  }

  // --- Phase 4: Finalize and Return ---

  report_stream << std::string(TOTAL_WIDTH, '=') << "\n";

  // Convert std::string to C-string (char*) for the C API
  std::string final_str = report_stream.str();
  char *cstr = new char[final_str.length() + 1];
  std::strcpy(cstr, final_str.c_str());
  return cstr;
}

/**
 * @brief Helper function to free the memory allocated by PDM_timer_get_report_string.
 * The C API requires the user to free memory allocated in the C++ layer.
 */
void PDM_timer_free_string(char* str) {
  delete[] str;
}


/**
 * @brief Prints the call tree to the console using the string builder function.
 */
void
PDM_timer_print
(
  PDM_timer_t *timer,
  int          mode
)
{
  char* report_str = PDM_timer_get_report_string(timer, mode);
  std::cout << report_str;
  PDM_timer_free_string(report_str); // Important: Free the dynamically allocated memory
}


/**
 * @brief Génère un dump JSON de l'arbre de calls (Local seulement).
 */
void
PDM_timer_dump_json
(
        PDM_timer_t *timer,
  const char        *filename
)
{
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_TIMER ERROR: Could not open file for JSON dump: %s \n", filename);
    return;
  }

  fprintf(fp, "{\n");
  fprintf(fp, "  \"profiling_data\": [\n");

  bool first_event = true;
  for (auto& pair : timer->root_event.children) {
    calculate_exclusive_and_dump(*pair.second, 1, fp, first_event);
  }

  fprintf(fp, "\n");
  fprintf(fp, "  ]\n");
  fprintf(fp, "}\n");

  fclose(fp);
}

/**
 * @brief Fonction récursive pour parcourir l'arbre et construire la map (Path -> Pointeur).
 * * @param node Le n?ud courant (commence à la racine ou à un enfant direct de la racine).
 * @param output_map La map à remplir.
 */
void collect_timer(_pdm_timer_event_t* node, std::map<std::string, _pdm_timer_event_t*>& output_map) {

  if (node->event_name != "__ROOT__") {
    output_map[node->path_name] = node;
  }

  // Parcours récursif de tous les enfants dans l'ordre d'insertion
  for (auto& child_name : node->child_insertion_order) {
    collect_timer(node->children.at(child_name).get(), output_map);
  }
}


int
get_serialize_timer_event_size
(
  _pdm_timer_event_t* node
)
{
  return (int) (sizeof(int) + 5 * sizeof(double));
}

void
PDM_timer_gather
(
  PDM_timer_t *timer
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(timer->comm, &i_rank);
  PDM_MPI_Comm_size(timer->comm, &n_rank);

  // Collect all data with tmp dict base on path
  std::map<std::string, _pdm_timer_event_t*> lflat_timer;
  collect_timer(&(timer->root_event), lflat_timer);

  int n_send        = lflat_timer.size();
  int n_send_path   = 0;
  for (auto& path_and_timer : lflat_timer) {
    std::cout << path_and_timer.first << std::endl;
    n_send_path += path_and_timer.first.size()+1;
  }


  std::vector<char>   send_buffer_path(    n_send_path);
  std::vector<int>    send_buffer_data(    n_send      );
  std::vector<double> send_buffer_time(6 * n_send      );

  int idx_write = 0;
  for (auto& path_and_timer : lflat_timer) {
    auto& path = path_and_timer.first;

    send_buffer_path.insert(send_buffer_path.end(), path.begin(), path.end());
    send_buffer_path.push_back('\0');

    auto& event = path_and_timer.second;
    send_buffer_data[  idx_write  ] = event->n_call;
    send_buffer_time[6*idx_write  ] = event->t1;
    send_buffer_time[6*idx_write+1] = event->t_run_inclusive;
    send_buffer_time[6*idx_write+2] = event->t_run_exclusive;
    send_buffer_time[6*idx_write+3] = event->t_children_sum;
    send_buffer_time[6*idx_write+4] = event->t_sync_entry;
    send_buffer_time[6*idx_write+5] = event->t_sync_exit;
    idx_write++;
  }

  int *gn_send_data = PDM_array_zeros_int(n_rank);
  int *gn_send_time = PDM_array_zeros_int(n_rank);
  PDM_MPI_Gather(&n_send      , 1, PDM_MPI_INT,
                  gn_send_data, 1, PDM_MPI_INT,
                  0,
                  timer->comm);

  int *gn_send_path_data = PDM_array_zeros_int(n_rank);;
  PDM_MPI_Gather(&n_send_path      , 1, PDM_MPI_INT,
                  gn_send_path_data, 1, PDM_MPI_INT,
                  0,
                  timer->comm);

  int *gn_send_path_data_idx = NULL;
  int *gn_send_data_idx      = NULL;
  int *gn_send_time_idx      = NULL;

  int n_g_data_recv      = 0;
  int n_g_data_path_recv = 0;
  if(i_rank == 0) {
    gn_send_path_data_idx = PDM_array_new_idx_from_sizes_int(gn_send_path_data, n_rank);
    gn_send_data_idx      = PDM_array_new_idx_from_sizes_int(gn_send_data     , n_rank);
    gn_send_time_idx      = PDM_array_new_idx_from_sizes_int(gn_send_data     , n_rank);
    n_g_data_recv      = gn_send_data_idx     [n_rank];
    n_g_data_path_recv = gn_send_path_data_idx[n_rank];

    for(int i = 0; i < n_rank; ++i) {
      gn_send_time[i] = gn_send_data[i] * 6;
    }
    for(int i = 0; i < n_rank+1; ++i) {
      gn_send_time_idx[i] *= 6;
    }

  }

  std::vector<char>   g_path(    n_g_data_path_recv);
  std::vector<int>    g_data(    n_g_data_recv     );
  std::vector<double> g_time(6 * n_g_data_recv     );

  PDM_MPI_Gatherv(send_buffer_path.data(), 1, PDM_MPI_CHAR,
                  g_path          .data(), gn_send_path_data, gn_send_path_data_idx, PDM_MPI_CHAR,
                  0,
                  timer->comm);

  PDM_MPI_Gatherv(send_buffer_data.data(), n_send, PDM_MPI_INT,
                  g_data          .data(), gn_send_data, gn_send_data_idx, PDM_MPI_INT,
                  0,
                  timer->comm);

  PDM_MPI_Gatherv(send_buffer_time.data(), 6 * n_send, PDM_MPI_DOUBLE,
                  g_time          .data(), gn_send_time, gn_send_time_idx, PDM_MPI_DOUBLE,
                  0,
                  timer->comm);

  if(i_rank == 0) {
    std::cout << std::string(g_path.data()) << std::endl;
  }

  printf("n_send_path = %i \n", n_send_path);

  PDM_log_trace_array_int   (g_data.data(), 1 * n_g_data_recv, "g_data ::");
  PDM_log_trace_array_double(g_time.data(), 6 * n_g_data_recv, "g_time ::");
  PDM_log_trace_array_double(send_buffer_time.data(), 6 * n_send, "send_buffer_time ::");

  // Create flat profile and reduce all data
  int idx_read_path = 0;
  int idx_read      = 0;
  std::map<std::string, _pdm_global_stat_t> gflat_timer;
  for(int t_rank = 0; t_rank < n_rank; ++t_rank) {

    for(int i = 0; i < gn_send_data[i]; ++i) {
      const char* path_start = g_path.data() + idx_read_path;
      std::string path(path_start);
      idx_read_path += path.length() + 1;

      long ln_call = g_data[idx_read];
      // double lt1 = g_time[6*idx_read];
      double lt_run_inclusive = g_time[6*idx_read+1];
      double lt_run_exclusive = g_time[6*idx_read+2];
      double lt_children_sum  = g_time[6*idx_read+3];
      double lt_sync_entry    = g_time[6*idx_read+4];
      double lt_sync_exit     = g_time[6*idx_read+5];

      _pdm_global_stat_t& g_record = gflat_timer[path];

      g_record.n_call           += ln_call;
      g_record.t_sum_inclusive  += lt_run_inclusive;
      g_record.t_sum_exclusive  += lt_run_exclusive;
      g_record.t_sum_sync_entry += lt_sync_entry;
      g_record.t_sum_sync_exit  += lt_sync_exit;

      // MIN/MAX (t_run_inclusive)
      if(lt_run_inclusive < g_record.t_min_run_inclusive) {
        g_record.t_min_run_inclusive = lt_run_inclusive;
        g_record.rank_min_inclusive  = t_rank;
      }
      if(lt_run_inclusive > g_record.t_max_run_inclusive) {
        g_record.t_max_run_inclusive = lt_run_inclusive;
        g_record.rank_max_inclusive  = t_rank;
      }

      // MIN/MAX (t_run_exclusive)
      if(lt_run_exclusive < g_record.t_min_run_exclusive) {
        g_record.t_min_run_exclusive = lt_run_exclusive;
        g_record.rank_min_exclusive  = t_rank;
      }
      if(lt_run_exclusive > g_record.t_max_run_exclusive) {
        g_record.t_max_run_exclusive = lt_run_exclusive;
        g_record.rank_max_exclusive  = t_rank;
      }

      // MIN/MAX (t_sync_entry)
      if(lt_sync_entry < g_record.t_min_sync_entry) {
        g_record.t_min_sync_entry = lt_sync_entry;
        g_record.rank_min_sync_entry = t_rank;
      }
      if(lt_sync_entry > g_record.t_max_sync_entry) {
        g_record.t_max_sync_entry = lt_sync_entry;
        g_record.rank_max_sync_entry = t_rank;
      }

      // MIN/MAX (t_sync_exit)
      if(lt_sync_exit < g_record.t_min_sync_exit) {
        g_record.t_min_sync_exit = lt_sync_exit;
        g_record.rank_min_sync_exit = t_rank;
      }
      if(lt_sync_exit > g_record.t_max_sync_exit) {
        g_record.t_max_sync_exit = lt_sync_exit;
        g_record.rank_max_sync_exit = t_rank;
      }

      idx_read++;
    }
  }

  // All data is computed, finalise mean
  for (auto& pair : gflat_timer) {
    _pdm_global_stat_t& g_record = pair.second;
    g_record.t_mean_run_inclusive = g_record.t_sum_inclusive  / n_rank;
    g_record.t_mean_run_exclusive = g_record.t_sum_exclusive  / n_rank;
    g_record.t_mean_sync_entry    = g_record.t_sum_sync_entry / n_rank;
    g_record.t_mean_sync_exit     = g_record.t_sum_sync_exit  / n_rank;
  }



  PDM_free(gn_send_path_data    );
  PDM_free(gn_send_data         );
  PDM_free(gn_send_path_data_idx);
  PDM_free(gn_send_data_idx     );

  free(gn_send_path_data);
}


void
PDM_timer_free2
(
  PDM_timer_t *timer
)
{
  delete timer;
}


// OLD

/*----------------------------------------------------------------------------
 * Creation d'un objet timer
 *
 * return
 *   timer
 *
 *----------------------------------------------------------------------------*/
PDM_timer_t*
PDM_timer_create
(
  void
)
{
  PDM_timer_t *timer;

  PDM_malloc(timer, 1, PDM_timer_t);
  PDM_timer_init(timer);

  return timer;
}

/*----------------------------------------------------------------------------
 * Debut la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_init(PDM_timer_t *timer)
{
#if defined (PDM_HAVE_GETRUSAGE)
  timer->t_cpu_u = 0.;
  timer->t_cpu_s = 0.;
#endif
  timer->t_cpu = 0.;
  timer->t_elapsed = 0.;
  timer->indic = 0;
}

/*----------------------------------------------------------------------------
 * Reprend la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_resume(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_reprise : \n"
            "La mesure d'une tranche est deja en cours\n");
    exit(EXIT_FAILURE);
  }
#if defined (PDM_HAVE_GETRUSAGE)
 {
   struct rusage  usage;

   if (getrusage(RUSAGE_SELF, &usage) == 0) {
     timer->t_cpu_u_debut = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.e-6;
     timer->t_cpu_s_debut = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1.e-6;
     timer->t_cpu_debut   = timer->t_cpu_u_debut + timer->t_cpu_s_debut;
   }
 }
#else
  timer->t_cpu_debut = clock();
#endif
  gettimeofday(&(timer->t_elaps_debut), NULL);
  timer->indic = 1;
}

/*----------------------------------------------------------------------------
 * Suspend la mesure du temps ecoule et incremente le temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_hang_on(PDM_timer_t *timer)
{

  if (!timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_suspend : \n"
            "La mesure de temps n'a pas ete declenchee par PDM_timer_reprise\n");
    exit(EXIT_FAILURE);
  }

  /* Recuperation du temps CPU et elaps courant */

  struct timeval t_elaps_fin;
  gettimeofday(&t_elaps_fin, NULL);

  /* Ajout de la tranche mesuree au temps cumule */

  long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
                         (timer->t_elaps_debut.tv_usec + 1000000 *
                          timer->t_elaps_debut.tv_sec);

  double tranche_elapsed_max = (double) tranche_elapsed;
  timer->t_elapsed += tranche_elapsed_max/1000000.;

#if defined (PDM_HAVE_GETRUSAGE)
 {
   struct rusage  usage;

   if (getrusage(RUSAGE_SELF, &usage) == 0) {
     timer->t_cpu_u += usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.e-6 - timer->t_cpu_u_debut;
     timer->t_cpu_s += usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1.e-6 - timer->t_cpu_s_debut;
     timer->t_cpu    = timer->t_cpu_u + timer->t_cpu_s;
   }
 }
#else
  clock_t t_cpu_fin = clock();
  double tranche_cpu = (double) (t_cpu_fin - timer->t_cpu_debut);
  double tranche_cpu_max = tranche_cpu;
  timer->t_cpu += tranche_cpu_max/CLOCKS_PER_SEC;
#endif
  timer->indic = 0;
}

/*----------------------------------------------------------------------------
 * Retourne le temps CPU en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_cpu : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_cpu\n");
    exit(EXIT_FAILURE);
  }
  return timer->t_cpu;
}

/*----------------------------------------------------------------------------
 * Retourne le temps CPU user en secondes (-1 si indisponible)
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu_user(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_cpu_user : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_cpu\n");
    exit(EXIT_FAILURE);
  }
#if defined (PDM_HAVE_GETRUSAGE)
  return timer->t_cpu_u;
#else
  return -1.;
#endif
}

/*----------------------------------------------------------------------------
 * Retourne le temps CPU systeme en secondes (-1 si indisponible)
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu_sys(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_cpu_user : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_cpu\n");
    exit(EXIT_FAILURE);
  }
#if defined (PDM_HAVE_GETRUSAGE)
  return timer->t_cpu_s;
#else
  return -1.;
#endif
}

/*----------------------------------------------------------------------------
 * Retourne le temps elaps en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_elapsed(PDM_timer_t *timer)
{
  if (timer->indic) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur PDM_timer_get_elapsed : \n"
            "Mesure d'une tranche en cours : faire appel a PDM_timer_suspend avant "
            "PDM_timer_get_elapsed\n");
    exit(EXIT_FAILURE);
  }
  return timer->t_elapsed;
}

/*----------------------------------------------------------------------------
 * Destruction d'un objet timer
 *
 * parameters :
 *   timer            <-- Timer
 *
 *----------------------------------------------------------------------------*/

void PDM_timer_free(PDM_timer_t *timer)
{
  PDM_free(timer);
}







#ifdef __cplusplus
}
#endif /* __cplusplus */
