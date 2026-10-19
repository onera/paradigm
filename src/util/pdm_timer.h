#ifndef __PDM_TIMER_H__
#define __PDM_TIMER_H__

/*----------------------------------------------------------------------------*/
#include "pdm.h"
#include "pdm_mpi.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure for measuring execution times
 *----------------------------------------------------------------------------*/

typedef struct _pdm_timer_t PDM_timer_t;

/*============================================================================
 * Public function interfaces
 *============================================================================*/


/**
 * \brief Initializes and creates a new \ref PDM_timer_t context.
 *
 * This function allocates and returns the main structure (\ref PDM_timer_t) used to record, organize, and aggregate event profiling data within a parallel application.
 *
 * The timer architecture is based on a **hierarchical timing system** that uses a
 * call stack to capture parent-child relationships between events. This structure
 * allows reports to be visualized in a clear, indented tree format.
 *
 * \par Key Profiling Concepts
 * The timer records several critical time metrics:
 * - **Inclusive Time:**
 * Measures the total elapsed time (Wall-Clock Time, typically from \ref PDM_MPI_Wtime)
 * between the start and end of an event. This **includes** the time spent in any
 * child functions (events) it calls.
 * - **Exclusive Time:**
 * Measures the **net time** spent inside the function itself. It is calculated by
 * subtracting the inclusive time of all direct child events from the event's
 * Inclusive Time. This is the key metric for identifying a function's direct computational cost.
 * - **Synchronization-Aware Timing :**
 * In addition to inclusive/exclusive times, the timer captures the time spent in
 * synchronization operations (e.g., MPI barriers, waits) immediately **before** and **after** the event's main measurement.
 * These metrics are essential for diagnosing **load imbalance** or **waiting latency** * in a parallel context.
 *
 * \param [in] comm The MPI communicator (\ref PDM_MPI_Comm) on which the timer operates.
 * This communicator will be used for the global aggregation of results.
 *
 * \return A pointer to the newly allocated timer structure (\ref PDM_timer_t*).
 * The caller is responsible for freeing this structure using \ref PDM_timer_free.
 *
 */
PDM_timer_t*
PDM_timer_create
(
  PDM_MPI_Comm comm
);


/**
 * \brief Marks the start of a timed event, pushing it onto the timer's call stack.
 *
 * This function initiates the measurement for a named event, establishing its place
 * in the hierarchical structure of the profiling data. It automatically sets up the
 * parent-child relationship with the currently running event (if one exists).
 *
 * \par Hierarchical Event Management
 * When called, \ref PDM_timer_start does the following:
 * 1. **Pushes the Event:** The new event defined by \p name is pushed onto the timer's internal call stack.
 * 2. **Sets Parent:** The event that was previously at the top of the stack (if any) is automatically designated as the **parent**.
 * 3. **Records Path:** The unique **path name** (e.g., "Main/FunctionA/LoopB") is derived and stored, which is crucial for global aggregation.
 * 4. **Starts Timing:** The Wall-Clock Time (\ref PDM_MPI_Wtime) is recorded to mark the start of the event's **Inclusive Time**.
 *
 * \param[in] timer The timer context structure (\ref PDM_timer_t*) created by \ref PDM_timer_create.
 * \param[in] name The simple, descriptive name of the event (e.g., "step 1"). This name is appended to the parent's path name.
 * \param[in] force_synchro An integer flag controlling synchronization behavior at the event start:
 * - If non-zero : The function executes an **MPI Barrier** * on the timer's communicator just before starting the event measurement.
 * The time spent waiting at the barrier is recorded.
 * This helps isolate the impact of load imbalance on the event's start time.
 * - If zero : No synchronization is performed, and timing starts immediately.
 */
void
PDM_timer_start
(
        PDM_timer_t *timer,
  const char        *name,
        int          force_synchro
);

/**
 * \brief Marks the end of a timed event, calculating its duration and removing it from the call stack.
 *
 * This function concludes the timing for the event identified by \p name, which must match
 * the event currently at the top of the timer's internal call stack. It calculates the
 * event's total Inclusive Time and finalizes all associated metrics before popping it from the stack.
 *
 * \par Event Finalization Process
 * When called, \ref PDM_timer_end performs the following crucial steps:
 * 1. **Records End Time:** The Wall-Clock Time (\ref PDM_MPI_Wtime) is recorded to mark the end of the event.
 * 2. **Calculates Duration:** The total **Inclusive Time**is computed by taking the difference between the end and start times.
 * 3. **Updates Parent:** The newly calculated Inclusive Time is subtracted from the parent event's potential Exclusive Time (by summing it to the parent's children time total).
 * 4. **Pops Event:** The finished event is popped from the call stack, restoring the parent event to the top of the stack.
 *
 * \param[in] timer The timer context structure (\ref PDM_timer_t) created by \ref PDM_timer_create.
 * \param[in] name The name of the event being terminated. This name *must* exactly match the name of the event at the top of the stack. A mismatch indicates a logic error in the application's timer calls.
 * \param[in] force_synchro An integer flag controlling synchronization behavior upon event completion:
 * - If non-zero : The function executes an **MPI Barrier** on the timer's communicator just *after* recording the end time.
 * The time spent waiting at the barrier is recorded.
 * This helps assess the impact of waiting for other ranks to finish the same event.
 * - If zero : No synchronization is performed.
 *
 */
void
PDM_timer_end
(
        PDM_timer_t *timer,
  const char        *name,
        int          force_synchro
);

/**
 * \brief Get the duration an number of calls of an event
 *
 * \param [in]  timer    PDM_timer_t instance
 * \param [in]  path     Path of the event
 * \param [out] duration Duration (in seconds) of the event
 *
 * \return Number of event calls
 */
long
PDM_timer_get
(
        PDM_timer_t *timer,
  const char        *path,
        double      *duration
);

/**
 * \brief Aggregates data and prints the profiling report in current process to the standard output.
 * \par Report Modes
 * The report format is controlled by the mode parameter:
 * - \ref PDM_TIMER_REPORT_HIERARCHICAL: The report preserves the call structure, indenting child events
 * under their parents. This is ideal for analyzing call flow and identifying top-level consumers.
 * - \ref PDM_TIMER_REPORT_FLAT: The report lists all unique events alphabetically by their full path name
 * (e.g., "Main/Kernel/Loop"), ignoring the call hierarchy for a concise overview of every measured function.
 *
 * \param timer [in] The timer context structure (\ref PDM_timer_t*)
 * \param mode  [in] Report mode (\ref PDM_timer_report_t)
 */
void
PDM_timer_print
(
  PDM_timer_t        *timer,
  PDM_timer_report_t  mode
);

/**
 * \brief Aggregates data and log the profiling report in current process to the paradigm logger
 * \par Report Modes
 * The report format is controlled by the mode parameter:
 * - \ref PDM_TIMER_REPORT_HIERARCHICAL: The report preserves the call structure, indenting child events
 *         under their parents. This is ideal for analyzing call flow and identifying top-level consumers.
 * - \ref PDM_TIMER_REPORT_FLAT: The report lists all unique events alphabetically by their full path name
 * (e.g., "Main/Kernel/Loop"), ignoring the call hierarchy for a concise overview of every measured function.
 *
 * \param timer [in] The timer context structure (\ref PDM_timer_t*)
 * \param mode  [in] Report mode (\ref PDM_timer_report_t)
 */
void
PDM_timer_log
(
  PDM_timer_t        *timer,
  PDM_timer_report_t  mode
);

/**
 * \brief Exports the **local**, per-process profiling data to a JSON file.
 *
 * This function serializes the full hierarchical tree of timed events recorded by
 * the calling process into a specified JSON file. Since no aggregation
 * is performed, this report is crucial for **analyzing the timing behavior of
 * an individual process** and identifying process-specific performance anomalies
 * or load imbalance sources.
 *
 * \par JSON Structure and Content
 * The exported JSON file maintains the **exact hierarchical structure** established
 * during the runtime (parent-child relationships). Each event object includes the
 * following per-process metrics:
 * - **name:** The simple name of the event (e.g., "LoopA").
 * - **n_call:** The number of times this event was executed locally.
 * - **t_inclusive:** The total Wall-Clock Time for the event.
 * - **t_exclusive:** The net time spent within the event (calculated during the dump process).
 * - **t_sync_entry:** Synchronization time recorded at the start.
 * - **t_sync_exit:** Synchronization time recorded at the end.
 * - **children:** An array containing the full JSON objects of child events.
 *
 * \note This function is typically called by **all ranks**, but the output file will
 * contain only the local data of the calling process. To distinguish reports,
 * the \p filename should usually include the rank ID (e.g., "profiling_0.json",
 * "profiling_1.json").
 *
 * \param [in] timer    The timer context structure (\ref PDM_timer_t*) containing the local event tree.
 * \param [in] filename The path and name of the JSON file to be created.
 *
 */
void
PDM_timer_dump_json
(
        PDM_timer_t *timer,
  const char        *filename
);


/**
 * \brief Aggregates profiling data from all processes and writes the consolidated report to a file on rank 0.
 *
 * This function initiates the collective process of gathering local statistics from all ranks
 * into a single comprehensive set of **Global Statistics** (Mean, Min[R], Max[R]) on the
 * root process (rank 0).
 *
 * \par Report Format
 * The resulting file contains a structured, typically **text-based** report. This function is generally used to generate the final human-readable
 * terminal or log file report.
 *
 * \note This function is called **collectively**, but the final file output only occurs on the root rank (rank 0) after aggregation is complete.
 *
 * \param[in] timer    The timer context structure (\ref PDM_timer_t*).
 * \param[in] filename The path and name of the file to which the global aggregated data will be written (e.g., "pdm_report.txt").
 *
 */
void
PDM_timer_gather_dump
(
        PDM_timer_t *timer,
  const char        *filename
);

/**
 * \brief Aggregates profiling data and exports the **Global Statistics** in a hierarchical JSON format on rank 0.
 *
 * This function first performs the global aggregation (calculating Mean, Min, Max, etc., across all ranks).
 * It then generates a structured JSON file that uses the original event call tree (\ref PDM_timer_t hierarchy)
 * but populates the timing fields with the **aggregated global results**.
 *
 * \par JSON Structure and Content
 * The output JSON is highly machine-readable, maintaining the parent-child relationships while
 * exposing all global metrics for each event:
 * - **path_name:** Unique identifier (e.g., "Main/Kernel").
 * - **n_call:** Total calls across all processes.
 * - **t_inclusive_mean/min/max:** Global statistics for Inclusive Time.
 * - **t_exclusive_mean:** Global mean for Exclusive Time.
 * - **r_inclusive_min/max:** Rank IDs where the minimum and maximum times occurred.
 *
 * \param [in] timer    The timer context structure (\ref PDM_timer_t*).
 * \param [in] filename The path and name of the JSON file to be created (e.g., "global_profiling.json").
 *
 */
void
PDM_timer_gather_dump_json
(
        PDM_timer_t *timer,
  const char        *filename
);

/**
 * \brief Releases all dynamically allocated resources associated with the PDM timer context.
 *
 * This function is the final cleanup step for the timer. It recursively frees all
 * memory used by the timer structure, including the entire event tree, the call stack,
 * and any globally aggregated statistics that may have been stored.
 *
 * \note This function should be called **after** all profiling and reporting is complete.
 * **Failure to call this function will result in memory leaks.**
 *
 * \param [in] timer The timer context structure (\ref PDM_timer_t*) to be destroyed and freed.
 */
void
PDM_timer_free
(
  PDM_timer_t *timer
);

/**
 * \brief Generates and returns the local, per-process profiling report as a dynamically allocated string.
 *
 * This utility function calls the internal reporting engine to format the local event
 * data into a human-readable text report. This report reflects the exact times and
 * call stack structure observed by the calling process (rank).
 *
 * \par Memory Management
 * The function returns a newly allocated C-style string (`char*`). **The caller is
 * responsible for freeing this memory** using the appropriate deallocation function.
 *
 * \par Report Modes
 * The structure of the generated report is controlled by the \p mode parameter:
 * - \ref PDM_TIMER_REPORT_HIERARCHICAL: The report displays the event timing data in a hierarchical,
 * indented format, mirroring the call stack of the process.
 * - \ref PDM_TIMER_REPORT_FLAT: The report lists all unique events alphabetically by their full path
 * name, providing a flat view of the function costs.
 *
 * \param [in] timer The timer context structure (\ref PDM_timer_t*) containing the local event tree.
 * \param [in] mode  Report mode (\ref PDM_timer_report_t)
 *
 * \return A dynamically allocated C-string containing the formatted local report. Returns NULL or an empty string on allocation failure.
 *
 */
char*
PDM_timer_get_report_string
(
  PDM_timer_t        *timer,
  PDM_timer_report_t  mode
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_TIMER_H__ */
