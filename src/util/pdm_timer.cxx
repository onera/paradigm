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

#include <unordered_map>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <functional>
#include <iomanip>

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

  // Order list: Stores keys in the order they were first inserted (crucial for Call Tree sequence)
  std::vector<std::string> child_insertion_order;

  // Tree Structure: std::unordered_map for O(1) lookup and ownership
  std::unordered_map<std::string, std::unique_ptr<_pdm_timer_event_t>> children;

  // Node Statistics
  long   n_call = 0;
  double t1     = 0.;             // Start time of the active call
  double t_run_inclusive = 0.;    // Total inclusive time
  double t_run_exclusive = 0.;    // Exclusive time (calculated at dump)
  double t_children_sum  = 0.;    // Sum of children's inclusive times

  // Recursion Management
  int    is_active_count = 0;

  _pdm_timer_event_t() = default;
  _pdm_timer_event_t(_pdm_timer_event_t const&) = delete;
  _pdm_timer_event_t& operator=(_pdm_timer_event_t const&) = delete;
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
 * @brief Traverses the tree for hierarchical and aligned console printing.
 */
void traverse_and_print_aligned(_pdm_timer_event_t& node, int depth, size_t name_width, int time_width, int ncall_width) {

  if (node.event_name == "__ROOT__") {
    for (auto& child_name : node.child_insertion_order) {
      traverse_and_print_aligned(*node.children.at(child_name), depth, name_width, time_width, ncall_width);
    }
    return;
  }

  // 1. Calculate Exclusive Time
  double children_time_sum = 0.0;
  for (auto& child_name : node.child_insertion_order) {
    children_time_sum += node.children.at(child_name)->t_run_inclusive;
  }
  double t_exclusive = node.t_run_inclusive - children_time_sum;

  // 2. Display

  // --- Column 1: Event Name (with Indentation) ---
  std::string indent = "";
  for (int i = 0; i < depth; ++i) {
    indent += "  |";
  }
  std::string indented_name = indent + node.event_name;

  std::cout << std::left << std::setw(name_width) << indented_name;

  // --- Column 2: N_Call ---
  std::cout << std::right << std::setw(ncall_width) << node.n_call;

  // --- Column 3: T_Inclusive ---
  std::cout << " |"
            << std::right << std::setw(time_width)
            << std::fixed << std::setprecision(6)
            << node.t_run_inclusive;

  // --- Column 4: T_Exclusive ---
  std::cout << " |"
            << std::right << std::setw(time_width)
            << std::fixed << std::setprecision(6)
            << t_exclusive;

  std::cout << std::endl;

  // 3. Recurse (using the insertion order)
  for (auto& child_name : node.child_insertion_order) {
    traverse_and_print_aligned(*node.children.at(child_name), depth + 1, name_width, time_width, ncall_width);
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
    // Iterate using the guaranteed insertion order
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
    fprintf(fp, "\"t_exclusive\": %12.5e", node.t_run_exclusive);

    // Add children if present
    if (!node.children.empty()) {
      fprintf(fp, ",\n");
      for (int i = 0; i < (indent * 2) + 2; ++i) fprintf(fp, " ");
      fprintf(fp, "\"children\": [\n");

      bool current_level_first_child = true;
      // Recursive call using the insertion order
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

  if (force_synchro == 1) {
    PDM_MPI_Barrier(timer->comm);
  }

  _pdm_timer_event_t* parent_node = get_active_parent(timer);
  _pdm_timer_event_t* event_ptr;

  // Check existence locally in the parent's children map.
  if (parent_node->children.count(current_name) == 0) {

    // --- CRUCIAL: Record the insertion order first ---
    parent_node->child_insertion_order.push_back(current_name);

    // Allocation via unique_ptr
    parent_node->children[current_name] = std::make_unique<_pdm_timer_event_t>();
    event_ptr = parent_node->children.at(current_name).get();

    // Initialization
    event_ptr->event_name = current_name;
    event_ptr->parent_name = parent_node->event_name;
  } else {
    // If the child exists (repeated or looped call), retrieve it.
    event_ptr = parent_node->children.at(current_name).get();
  }

  _pdm_timer_event_t& event = *event_ptr;

  // Time and Stack Management
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

  // 2. Retirer l'événement de la pile
  timer->call_stack.pop_back();

  // 3. Récupérer le n?ud parent (le nouveau sommet de la pile)
  _pdm_timer_event_t* parent_node = get_active_parent(timer);

  if (force_synchro == 1) {
    PDM_MPI_Barrier(timer->comm);
  }
  current_node->is_active_count--;

  if (current_node->is_active_count == 0) { // Fin du bloc externe
    double dt = PDM_MPI_Wtime() - current_node->t1;
    current_node->t_run_inclusive += dt;

    // 4. Mettre à jour le temps des enfants pour le parent (nécessaire pour le calcul exclusif du parent)
    parent_node->t_children_sum += dt;
  }
}


void PDM_timer_print(PDM_timer_t *timer, int mode) {
  if (timer->root_event.children.empty()) {
    std::cout << "PDM_TIMER: No events recorded." << std::endl;
    return;
  }

  // --- Phase 1: Calculate Dynamic Widths ---
  size_t max_name_width = 11; // Min width for "Event Name"

  if (mode == 0) { // Hierarchical Mode
    calculate_max_widths(timer->root_event, 0, max_name_width);
  } else { // Flat/Raw Mode
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

  // Fixed widths for numeric columns
  const int N_CALL_COL_WIDTH = 10;
  const int TIME_COL_WIDTH = 18;

  const int TOTAL_WIDTH = max_name_width + N_CALL_COL_WIDTH + (TIME_COL_WIDTH * 2) + 4;

  // --- Phase 2: Print Header and Content ---

  std::cout << "\n" << std::string(TOTAL_WIDTH, '=') << "\n";
  std::cout << "PDM TIMER REPORT (" << (mode == 1 ? "FLAT/RAW MODE" : "HIERARCHICAL MODE") << ")\n";
  std::cout << std::string(TOTAL_WIDTH, '=') << "\n";

  // Column Header
  std::cout << std::left << std::setw(max_name_width) << "Event Name (Path)";
  std::cout << std::right << std::setw(N_CALL_COL_WIDTH) << "Calls";
  std::cout << " |" << std::right << std::setw(TIME_COL_WIDTH) << "Time Inclusive (s)";
  std::cout << " |" << std::right << std::setw(TIME_COL_WIDTH) << "Time Exclusive (s)";
  std::cout << "\n" << std::string(TOTAL_WIDTH, '-') << "\n";

  if (mode == 0) { // Hierarchical Mode
    for (auto& child_name : timer->root_event.child_insertion_order) {
      traverse_and_print_aligned(*timer->root_event.children.at(child_name), 0, max_name_width, TIME_COL_WIDTH, N_CALL_COL_WIDTH);
    }
  } else { // Flat/Raw Mode
    std::vector<_pdm_timer_event_t*> all_nodes;
    collect_all_nodes(&timer->root_event, all_nodes);

    for(auto* node : all_nodes) {
      // Recalculate stats
      double children_time_sum = 0.0;
      for (auto& child_name : node->child_insertion_order) {
        children_time_sum += node->children.at(child_name)->t_run_inclusive;
      }
      double t_exclusive = node->t_run_inclusive - children_time_sum;

      // Display
      std::cout << std::left << std::setw(max_name_width) << node->event_name;
      std::cout << std::right << std::setw(N_CALL_COL_WIDTH) << node->n_call;
      std::cout << " |"
                << std::right << std::setw(TIME_COL_WIDTH)
                << std::fixed << std::setprecision(6)
                << node->t_run_inclusive;
      std::cout << " |"
                << std::right << std::setw(TIME_COL_WIDTH)
                << std::fixed << std::setprecision(6)
                << t_exclusive;
      std::cout << std::endl;
    }
  }

  std::cout << std::string(TOTAL_WIDTH, '=') << "\n";
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
