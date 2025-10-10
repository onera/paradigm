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
  int exclusive_time_computed = 0;
  _pdm_timer_event_t root_event;
  std::vector<_pdm_timer_event_t*> call_stack;

  int is_gather = 0;
  std::map<std::string, _pdm_global_stat_t> gflat_timer;

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

static
void
collect_timer
(
  _pdm_timer_event_t*                         node,
  std::map<std::string, _pdm_timer_event_t*>& output_map
)
{

  if (node->event_name != "__ROOT__") {
    output_map[node->path_name] = node;
  }

  // Parcours récursif de tous les enfants dans l'ordre d'insertion
  for (auto& child_name : node->child_insertion_order) {
    collect_timer(node->children.at(child_name).get(), output_map);
  }
}

static
std::string
format_condensed_value
(
  double time,
  int    rank
)
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(4) << time;
  if (rank != -1) {
    ss << "[" << rank << "]";
  }
  return ss.str();
}

static
std::string
format_full_condensed_stat
(
  double t_mean,
  double t_min,
  int    r_min,
  double t_max,
  int    r_max,
  int    width
)
{
  std::string mean_str = format_condensed_value(t_mean, -1);
  std::string min_str  = format_condensed_value(t_min, r_min);
  std::string max_str  = format_condensed_value(t_max, r_max);

  std::string content = mean_str + "/" + min_str + "/" + max_str;

  // Assurer l'alignement à droite dans la largeur donnée
  std::stringstream ss_final;
  ss_final << std::right << std::setw(width) << content;
  return ss_final.str();
}

std::string
format_timer_line_main
(
       _pdm_timer_event_t&                        node,
 const std::string&                               indented_name,
       size_t                                     name_width,
       int                                        time_width,
       int                                        ncall_width,
 const std::map<std::string, _pdm_global_stat_t>* global_stats
)
{
  std::stringstream ss;
  const _pdm_global_stat_t* g_rec = nullptr;
  bool is_global_report = false;

  if (global_stats && global_stats->count(node.path_name)) {
    g_rec = &global_stats->at(node.path_name);
    is_global_report = true;
  }

  // --- 1. Event Name & Calls ---
  ss << std::left << std::setw(name_width) << indented_name;

  long n_call_to_display = is_global_report ? g_rec->n_call : node.n_call;
  ss << std::right << std::setw(ncall_width) << n_call_to_display;

  ss << " |"; // Séparateur

  // --- 2. Basculement Local vs Global ---
  if (is_global_report) {

    // --- COLONNES GLOBAL COMPACTES (MEAN/MIN[R]/MAX[R]) ---
    // T_Inclusive
    ss << std::right << format_full_condensed_stat(g_rec->t_mean_run_inclusive,
                                                   g_rec->t_min_run_inclusive, g_rec->rank_min_inclusive,
                                                   g_rec->t_max_run_inclusive, g_rec->rank_max_inclusive,
                                                   time_width);

    // T_Exclusive
    ss << " |" << std::right << format_full_condensed_stat(g_rec->t_mean_run_exclusive,
                                                           g_rec->t_min_run_exclusive, g_rec->rank_min_exclusive,
                                                           g_rec->t_max_run_exclusive, g_rec->rank_max_exclusive,
                                                           time_width);

    // T_Sync_Entry
    ss << " |" << std::right << format_full_condensed_stat(g_rec->t_mean_sync_entry,
                                                           g_rec->t_min_sync_entry, g_rec->rank_min_sync_entry,
                                                           g_rec->t_max_sync_entry, g_rec->rank_max_sync_entry,
                                                           time_width);

    // T_Sync_Exit
    ss << " |" << std::right << format_full_condensed_stat(g_rec->t_mean_sync_exit,
                                                           g_rec->t_min_sync_exit, g_rec->rank_min_sync_exit,
                                                           g_rec->t_max_sync_exit, g_rec->rank_max_sync_exit,
                                                           time_width);

  } else {

    // --- COLONNES LOCALES (Standard - 4 colonnes) ---
    double t_exclusive_local = node.t_run_inclusive - node.t_children_sum;
    if (t_exclusive_local < 0) t_exclusive_local = 0;

    // T_Inclusive
    ss << std::right << std::setw(time_width) << std::fixed << std::setprecision(6) << node.t_run_inclusive;

    // T_Exclusive
    ss << " |" << std::right << std::setw(time_width) << std::fixed << std::setprecision(6) << t_exclusive_local;

    // T_Sync_Entry
    ss << " |" << std::right << std::setw(time_width) << std::fixed << std::setprecision(6) << node.t_sync_entry;

    // T_Sync_Exit
    ss << " |" << std::right << std::setw(time_width) << std::fixed << std::setprecision(6) << node.t_sync_exit;
  }

  return ss.str();
}

static
void
collect_all_nodes
(
  _pdm_timer_event_t*               node,
  std::vector<_pdm_timer_event_t*>& all_nodes
)
{
  if (node->event_name != "__ROOT__") {
    all_nodes.push_back(node);
  }
  for (auto& child_name : node->child_insertion_order) {
    collect_all_nodes(node->children.at(child_name).get(), all_nodes);
  }
}

static
void
traverse_and_add_lines
(
        _pdm_timer_event_t&                        node,
        int                                        depth,
        size_t                                     name_width,
        int                                        time_width,
        int                                        ncall_width,
  const std::map<std::string, _pdm_global_stat_t>* global_stats,
        std::vector<std::string>&                  lines
)
{
  // Ignorer le n?ud racine, mais continuer à parcourir ses enfants
  if (node.event_name == "__ROOT__") {
    for (auto& child_name : node.child_insertion_order) {
      traverse_and_add_lines(*node.children.at(child_name), depth, name_width, time_width, ncall_width, global_stats, lines);
    }
    return;
  }

  // 1. Préparer l'indentation
  std::string indent = "";
  for (int i = 0; i < depth; ++i) {
    indent += "  |";
  }
  std::string indented_name = indent + node.event_name;

  // 2. Formater et ajouter la seule ligne (Local ou Global Condensé: MEAN/MIN[R]/MAX[R])
  lines.push_back(
    format_timer_line_main(node, indented_name, name_width, time_width, ncall_width, global_stats)
  );

  // 3. Parcourir les enfants
  for (auto& child_name : node.child_insertion_order) {
    traverse_and_add_lines(*node.children.at(child_name), depth + 1, name_width, time_width, ncall_width, global_stats, lines);
  }
}


static
void
_recursive_compute_exclusive_time
(
  _pdm_timer_event_t& node
)
{
  if (node.event_name != "__ROOT__") {
    double children_time_sum = 0.0;
    for (auto& child_name : node.child_insertion_order) {
      children_time_sum += node.children.at(child_name)->t_run_inclusive;
    }
    node.t_run_exclusive = node.t_run_inclusive - children_time_sum;

    if (!node.children.empty()) {
      for (auto& child_name : node.child_insertion_order) {
        _recursive_compute_exclusive_time(*node.children.at(child_name));
      }
    }
  }
}

static
void
_compute_exclusive_time
(
  PDM_timer_t* timer
)
{
  if(timer->exclusive_time_computed == 0) {
    for (auto& pair : timer->root_event.children) {
      _recursive_compute_exclusive_time(*pair.second);
    }
    timer->exclusive_time_computed = 1;
  }
}


static
void
_timer_gather
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
    n_send_path += path_and_timer.first.size()+1;
  }

  std::vector<char> send_buffer_path;
  send_buffer_path.reserve(n_send_path);
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

  int *gn_send_path_data = PDM_array_zeros_int(n_rank);
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

  PDM_MPI_Gatherv(send_buffer_path.data(), n_send_path, PDM_MPI_CHAR,
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

  if(0 == 1) {
    printf("n_send_path = %i \n", n_send_path);
    PDM_log_trace_array_int   (g_data.data()          , 1 * n_g_data_recv, "g_data           ::");
    PDM_log_trace_array_double(g_time.data()          , 6 * n_g_data_recv, "g_time           ::");
    PDM_log_trace_array_double(send_buffer_time.data(), 6 * n_send       , "send_buffer_time ::");
  }

  // Create flat profile and reduce all data
  int idx_read_path = 0;
  int idx_read      = 0;
  std::map<std::string, _pdm_global_stat_t> gflat_timer;
  for(int t_rank = 0; t_rank < n_rank; ++t_rank) {

    for(int i = 0; i < gn_send_data[t_rank]; ++i) {
      const char* path_start = g_path.data() + idx_read_path;
      std::string path(path_start);
      idx_read_path += path.length() + 1;

      long ln_call = g_data[idx_read];
      // double lt1 = g_time[6*idx_read];
      double lt_run_inclusive = g_time[6*idx_read+1];
      double lt_run_exclusive = g_time[6*idx_read+2];
      // double lt_children_sum  = g_time[6*idx_read+3];
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
    // std::cout << "gflat_timer = " << pair.first << " -> " << pair.first.size() << std::endl;
    _pdm_global_stat_t& g_record = pair.second;
    g_record.t_mean_run_inclusive = g_record.t_sum_inclusive  / n_rank;
    g_record.t_mean_run_exclusive = g_record.t_sum_exclusive  / n_rank;
    g_record.t_mean_sync_entry    = g_record.t_sum_sync_entry / n_rank;
    g_record.t_mean_sync_exit     = g_record.t_sum_sync_exit  / n_rank;
  }

  PDM_free(gn_send_time);
  PDM_free(gn_send_path_data    );
  PDM_free(gn_send_data         );
  PDM_free(gn_send_path_data_idx);
  PDM_free(gn_send_data_idx     );
  PDM_free(gn_send_time_idx     );

  free(gn_send_path_data);
  timer->is_gather = 1;
  timer->gflat_timer = std::move(gflat_timer);
}

static
void
_dump_json
(
  _pdm_timer_event_t& node,
  int                 indent,
  FILE*               fp,
  bool&               first_child
)
{
  if (node.event_name != "__ROOT__") {
    if (!first_child) {
      fprintf(fp, ",\n");
    }
    first_child = false;

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
        _dump_json(*node.children.at(child_name), indent + 2, fp, current_level_first_child);
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


static
std::stringstream
_pdm_timer_generate_report
(
        PDM_timer_t*                               timer,
        int                                        mode,
  const std::map<std::string, _pdm_global_stat_t>* global_stats,
  const std::string&                               report_title
)
{
  std::stringstream report_stream;

  // --- Phase 1: Détermination du contexte et des largeurs ---
  bool is_global_report = (global_stats != nullptr);

  size_t max_name_width = 11;
  calculate_max_widths(timer->root_event, 0, max_name_width);
  max_name_width = std::min(max_name_width + 2, (size_t)80);

  const int N_CALL_COL_WIDTH = 10;
  const int TIME_COL_WIDTH   = 35;
  const int NUM_TIME_COLS    = 4; // T_INC, T_EXC, T_SE, T_SX

  const int TOTAL_WIDTH = max_name_width + N_CALL_COL_WIDTH + (TIME_COL_WIDTH * NUM_TIME_COLS) + (NUM_TIME_COLS * 2) + 2;

  // --- Phase 2: Construction de l'En-tête ---
  report_stream << "\n" << std::string(TOTAL_WIDTH, '=') << "\n";
  report_stream << "PDM TIMER REPORT (" << report_title << ", "
                << (mode == 1 ? "FLAT/RAW MODE" : "HIERARCHICAL MODE") << ")\n";
  report_stream << std::string(TOTAL_WIDTH, '=') << "\n";

  // En-tête des colonnes
  report_stream << std::left << std::setw(max_name_width) << "Event Name (Path)";
  report_stream << std::right << std::setw(N_CALL_COL_WIDTH) << "Calls";

  // En-têtes des 4 métriques de temps
  if (is_global_report) {
    // En-tête condensé
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_INCLUSIVE  MEAN/MIN[R]/MAX[R] (s)";
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_EXCLUSIVE  MEAN/MIN[R]/MAX[R] (s)";
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_SYNC_ENTRY MEAN/MIN[R]/MAX[R] (s)";
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_SYNC_EXIT  MEAN/MIN[R]/MAX[R] (s)";
  } else {
    // En-tête local standard
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_INCLUSIVE  (s)";
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_EXCLUSIVE  (s)";
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_SYNC_ENTRY (s)";
    report_stream << " |" << std::right << std::setw(TIME_COL_WIDTH) << "T_SYNC_EXIT  (s)";
  }

  report_stream << "\n" << std::string(TOTAL_WIDTH, '-') << "\n";

  // --- Phase 3: Parcours et Ajout des Lignes ---
  std::vector<std::string> content_lines;

  if (mode == 0) { // Hierarchical Mode
    for (auto& child_name : timer->root_event.child_insertion_order) {
      traverse_and_add_lines(*timer->root_event.children.at(child_name), 0, max_name_width, TIME_COL_WIDTH, N_CALL_COL_WIDTH, global_stats, content_lines);
    }
  } else { // Flat/Raw Mode
    std::vector<_pdm_timer_event_t*> all_nodes;
    collect_all_nodes(&timer->root_event, all_nodes);

    // La largeur de l'espace vide pour aligner les stats MIN/MAX en mode plat
    const int EMPTY_COL_WIDTH = max_name_width + N_CALL_COL_WIDTH + 2;
    std::string empty_prefix = std::string(EMPTY_COL_WIDTH, ' ');

    for(auto* node : all_nodes) {
      content_lines.push_back(format_timer_line_main(*node,
                                                     node->event_name,
                                                     max_name_width,
                                                     TIME_COL_WIDTH,
                                                     N_CALL_COL_WIDTH,
                                                     global_stats));
    }
  }

  for(const auto& line : content_lines) {
    report_stream << line << "\n";
  }

  // --- Phase 4: Finalisation ---
  report_stream << std::string(TOTAL_WIDTH, '=') << "\n";

  return report_stream;
}

static
void
dump_global_hierarchical_json
(
       _pdm_timer_event_t&                        node,
       int                                        indent,
       FILE*                                      fp,
       bool&                                      first_child,
 const std::map<std::string, _pdm_global_stat_t>* global_stats // Pointeur vers les stats agrégées
)
{
  if (node.event_name != "__ROOT__") {
    if (!first_child) {
      fprintf(fp, ",\n");
    }
    first_child = false;

    // Trouver les statistiques globales correspondantes (g_rec)
    const _pdm_global_stat_t* g_rec = nullptr;
    bool is_global_found = false;

    if (global_stats) {
      auto it = global_stats->find(node.path_name);
      if (it != global_stats->end()) {
        g_rec = &(it->second);
        is_global_found = true;
      }
    }

    // --- JSON Serialization ---
    std::string current_indent(indent * 2, ' ');
    std::string inner_indent((indent * 2) + 2, ' ');

    fprintf(fp, "%s{\n", current_indent.c_str());

    // Nom
    fprintf(fp, "%s\"name\": \"%s\",\n", inner_indent.c_str(), node.event_name.c_str());
    fprintf(fp, "%s\"path_name\": \"%s\",\n", inner_indent.c_str(), node.path_name.c_str());

    // Utilisation des données agrégées (g_rec)
    long n_call = is_global_found ? g_rec->n_call : node.n_call;
    fprintf(fp, "%s\"n_call\": %ld,\n", inner_indent.c_str(), n_call);

    // Exportation des 14 valeurs globales (Mean/Min/Max pour 4 métriques)

    // T_Inclusive
    fprintf(fp, "%s\"t_inclusive_mean\": %12.5e,\n" , inner_indent.c_str(), is_global_found ? g_rec->t_mean_run_inclusive : 0.0);
    fprintf(fp, "%s\"t_inclusive_min\": %12.5e,\n"  , inner_indent.c_str(), is_global_found ? g_rec->t_min_run_inclusive  : 0.0);
    fprintf(fp, "%s\"r_inclusive_min\": %d,\n"      , inner_indent.c_str(), is_global_found ? g_rec->rank_min_inclusive   : -1 );
    fprintf(fp, "%s\"t_inclusive_max\": %12.5e,\n"  , inner_indent.c_str(), is_global_found ? g_rec->t_max_run_inclusive  : 0.0);
    fprintf(fp, "%s\"r_inclusive_max\": %d,\n"      , inner_indent.c_str(), is_global_found ? g_rec->rank_max_inclusive   : -1 );
    fprintf(fp, "%s\"t_exclusive_mean\": %12.5e,\n" , inner_indent.c_str(), is_global_found ? g_rec->t_mean_run_exclusive : 0.0);
    fprintf(fp, "%s\"t_sync_entry_mean\": %12.5e,\n", inner_indent.c_str(), is_global_found ? g_rec->t_mean_sync_entry    : 0.0);
    fprintf(fp, "%s\"t_sync_exit_mean\": %12.5e"    , inner_indent.c_str(), is_global_found ? g_rec->t_mean_sync_exit     : 0.0);

    // Add children if present
    if (!node.children.empty()) {
      // Rajouter une virgule si on a des enfants
      fprintf(fp, ",\n");
      fprintf(fp, "%s\"children\": [\n", inner_indent.c_str());

      bool current_level_first_child = true;
      for (auto& child_name : node.child_insertion_order) {
        dump_global_hierarchical_json(*node.children.at(child_name), indent + 1, fp, current_level_first_child, global_stats);
      }
      fprintf(fp, "\n");
      fprintf(fp, "%s]", inner_indent.c_str());
    }

    fprintf(fp, "\n");
    fprintf(fp, "%s}", current_indent.c_str());
  }
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

PDM_timer_t*
PDM_timer_create
(
  PDM_MPI_Comm comm
)
{
  PDM_timer_t *timer = new PDM_timer_t(); // This is done like this to emulate C and have portability with C

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


char*
PDM_timer_get_report_string
(
  PDM_timer_t *timer,
  int          mode
)
{
  int i_rank = -1;
  PDM_MPI_Comm_rank(timer->comm, &i_rank);

  std::string title = "LOCAL Rank " + std::to_string(i_rank);
  std::stringstream report_stream = _pdm_timer_generate_report(timer, mode, nullptr, title);
  std::string       final_str     = report_stream.str();

  char *cstr = NULL;
  PDM_malloc(cstr, final_str.length() + 1, char);
  strcpy(cstr, final_str.c_str());
  return cstr;
}


void
PDM_timer_print
(
  PDM_timer_t *timer,
  int          mode
)
{
  _compute_exclusive_time(timer);
  char* report_str = PDM_timer_get_report_string(timer, mode);
  std::cout << report_str;
  PDM_free(report_str);
}


void
PDM_timer_log
(
  PDM_timer_t *timer,
  int          mode
)
{
  _compute_exclusive_time(timer);
  char* report_str = PDM_timer_get_report_string(timer, mode);
  log_trace("%s", report_str);
  PDM_free(report_str);
}


void
PDM_timer_dump_json
(
        PDM_timer_t *timer,
  const char        *filename
)
{
  _compute_exclusive_time(timer);

  FILE *fp = fopen(filename, "w");
  if (!fp) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_TIMER ERROR: Could not open file for JSON dump: %s \n", filename);
    return;
  }

  fprintf(fp, "{\n");
  fprintf(fp, "  \"profiling_data\": [\n");

  bool first_event = true;
  for (auto& pair : timer->root_event.children) {
    _dump_json(*pair.second, 1, fp, first_event);
  }

  fprintf(fp, "\n");
  fprintf(fp, "  ]\n");
  fprintf(fp, "}\n");

  fclose(fp);
}


void
PDM_timer_gather_dump_json
(
        PDM_timer_t *timer,
  const char        *filename
)
{
  if(timer->is_gather == 0) {
    _compute_exclusive_time(timer);
    _timer_gather(timer);
  }

  int i_rank = -1;
  PDM_MPI_Comm_rank(timer->comm, &i_rank);

  if(i_rank != 0) {
    return; // Only rank 0 export data
  }

  // A completer ici
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "PDM_TIMER ERROR: Could not open file for JSON dump: %s \n", filename);
    return;
  }

  // 2. Écrire l'en-tête JSON
  fprintf(fp, "{\n");
  fprintf(fp, "  \"profiling_data_aggregated\": [\n");

  // 3. Lancer la récursion sur l'arbre local avec les stats globales
  bool first_event = true;
  for (auto& child_name : timer->root_event.child_insertion_order) {
    dump_global_hierarchical_json(*timer->root_event.children.at(child_name), 1, fp, first_event, &timer->gflat_timer);
  }

  // 4. Écrire le pied de page JSON
  fprintf(fp, "\n  ]\n");
  fprintf(fp, "}\n");

  fclose(fp);
}


void
PDM_timer_gather_dump
(
  PDM_timer_t *timer,
  char        *filename
)
{
  if(timer->is_gather == 0) {
    _compute_exclusive_time(timer);
    _timer_gather(timer);
  }

  int i_rank = -1;
  PDM_MPI_Comm_rank(timer->comm, &i_rank);
  if(i_rank == 0) {
    std::string title = "AGGREGATED GLOBAL (All Ranks)";
    std::stringstream report_stream = _pdm_timer_generate_report(timer, 0, &timer->gflat_timer, title);
    std::string       final_str     = report_stream.str();

    FILE *fp = fopen(filename, "w");
    if (!fp) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_TIMER ERROR: Could not open file for dump: %s \n", filename);
      return;
    }
    fprintf(fp, "%s", final_str.c_str());
    fclose(fp);
  }
}


void
PDM_timer_free
(
  PDM_timer_t *timer
)
{
  delete timer;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
