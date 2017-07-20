#ifndef PARALLELISM_HPP
#define PARALLELISM_HPP

#include <pthread.h>
#include <list>
#include <iostream>

#include "distributions.hpp" 
// todo: think of a better way to have multiple Rng's

extern pthread_mutex_t coutMutex;

/* todo: 
 * (1) when jobs are easy, it might be better to pass multiple
 * jobs at once (i.e. in serial) to the workers. Maybe this 
 * could be determined automatically?
 * (2) add a method addJobs that accepts ranges.
 * (3) re-implement this with C++11, or boost.
 */

/* the interface Job is used by WorkerPool, which contains 
 * a list of Jobs. Classes that inherit Job, can be handled
 * by WorkerPool, and must implement the method "execute"
 */
class Job {
	/* WorkerPool::workerThreadEntry needs to call setFinished()
	 * after Job::execute returns. But, I don't want anybody else to
	 * call setFinished! todo: better method?
	 */
	friend class WorkerPool;
public:
	Job(); // init private members
	virtual ~Job(); // obligatory virtual destructor
	bool isFinished();
	void waitForJobToFinish();
	void setDeleteWhenFinished();
protected:
	/* WorkerPool::workerThreadEntry calls this method.
	 * Maybe replace the Rng with a void* for portability?
	 */
	virtual bool execute(Rng & )=0;
private: // children don't need to mess with these members...
	bool finished;
	pthread_mutex_t finishedMutex;
	pthread_cond_t finishedConditionalVar;
	void setFinished();
	int priority; // todo..
	bool deleteWhenFinished;
};

class WorkerPool {
public:
	/* constuctor */
	WorkerPool();
	/* destructor */
	virtual ~WorkerPool();
	/* Initiate private members and set the number of workers.
	 * Returns true if the threads were successfully started,
	 * false if there was an error starting a thread
	 */
	bool initWorkerPool(int , Rng & );
	/* Will not return until the internal threads have exited. */
	void waitForWorkerPoolToExit();
	/* Wait for all Workers to finish their jobs. The Workers
	 * will then be waiting for new jobs. (So no joining...)
	 */
	void syncWorkerThreads();
	/* Send the worker pools a signal that there will be no more jobs.
	 * This will make them exit their WorkerThreadEntry function.
	 */
	void sendNoMoreJobsSignal();
	/* add a new job to the list of jobs.
	 */
	void addNewJob(Job* );
protected:
	// no protected members?
private:
	void workerThreadEntry(unsigned long ); // pass a seed
	static void* workerThreadEntryFunc(void* ); // the argument is WorkerArg* (c/c++ hack)
	pthread_mutex_t jobsMutex;
	pthread_cond_t jobsConditionalVar;
	pthread_mutex_t jobsyncMutex;
	pthread_cond_t jobsyncConditionalVar;
	volatile int jobsyncCounter; // todo: when to use the keyword volatile?
	std::list<Job*> jobs;
	volatile bool noMoreJobsFlag;
	pthread_t* workerThreads;
	pthread_attr_t attr;
	int numOfWorkerThreads;
	/* Upon creation of a thread, an argument can be passed. This argument
	 * needs to contain a pointer to this object (WorkerPoolWrapper),
	 * but we can add other stuff... In this case each worker gets a 
	 * seed for a thread-local rng.
	 */ 
	struct WorkerArg {
		// inline constructors
		WorkerArg() : wp(NULL), seed(0) {}
		WorkerArg(WorkerPool* wp, unsigned long seed) : wp(wp), seed(seed) {}
		// the actual data...
		WorkerPool* wp;
		unsigned long seed;
	};
	WorkerArg* workerArgs;
};

#endif
