#include "parallelism.hpp"

pthread_mutex_t coutMutex = PTHREAD_MUTEX_INITIALIZER;

// methods for Job

Job::Job() {
	finished = false;
	pthread_mutex_init(&finishedMutex, NULL);
	pthread_cond_init(&finishedConditionalVar, NULL);
	priority = 0; // todo: implement
	deleteWhenFinished = false;
}
Job::~Job() { 
	pthread_mutex_destroy(&finishedMutex);
	pthread_cond_destroy(&finishedConditionalVar);
}
bool Job::isFinished() {
	bool answer;
	pthread_mutex_lock(&finishedMutex);
	answer = finished;
	pthread_mutex_unlock(&finishedMutex);
	// notice that finished could have switched value before answer is returned...
	return answer; 
}
void Job::waitForJobToFinish() {
	pthread_mutex_lock(&finishedMutex);
	while ( !finished ) {
		pthread_cond_wait(&finishedConditionalVar, &finishedMutex);
	}
	pthread_mutex_unlock(&finishedMutex);
}
void Job::setFinished() {
	pthread_mutex_lock(&finishedMutex);
	finished = true;
	pthread_cond_broadcast(&finishedConditionalVar); // wake threads that are waiting for job to finish
	pthread_mutex_unlock(&finishedMutex);
}
void Job::setDeleteWhenFinished() {
	deleteWhenFinished = true;
}


// methods for WorkerPool

WorkerPool::WorkerPool() {
	numOfWorkerThreads = 0;
	noMoreJobsFlag = true;
	jobsyncCounter = 0;
	workerThreads = NULL;
	workerArgs = NULL;
	// init pthread goodies
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_mutex_init(&jobsMutex, NULL);
	pthread_mutex_init(&jobsyncMutex, NULL);
	pthread_cond_init(&jobsConditionalVar, NULL);
	pthread_cond_init(&jobsyncConditionalVar, NULL);
}
WorkerPool::~WorkerPool() {
	// assumes that someone has joined the threads. todo: make fool proof.
	delete[] workerThreads;
	delete[] workerArgs;
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&jobsMutex);
	pthread_mutex_destroy(&jobsyncMutex);
	pthread_cond_destroy(&jobsConditionalVar);
	pthread_cond_destroy(&jobsyncConditionalVar);
}
bool WorkerPool::initWorkerPool(int numOfWorkerThreads, Rng & rng) {
	this->numOfWorkerThreads = numOfWorkerThreads;
	workerThreads = new pthread_t[numOfWorkerThreads];
	noMoreJobsFlag = false;
	jobsyncCounter = 0;
    
    workerArgs = new WorkerArg[numOfWorkerThreads];
    for ( int j = 0; j < numOfWorkerThreads; ++j ) {
		workerArgs[j] = WorkerArg(this, rng.integer());
	}
	
   	if ( numOfWorkerThreads < 1 ) return false;
	else {
		bool ok = true;
		for ( int i = 0; i < numOfWorkerThreads; i++ ) {
			ok = ok && ( pthread_create(workerThreads+i, &attr, 
			                workerThreadEntryFunc, workerArgs+i) == 0 );
		}
		return ok;
	}
}
void WorkerPool::waitForWorkerPoolToExit() {
	sendNoMoreJobsSignal(); // otherwise the workers will not return
	for ( int i = 0; i < numOfWorkerThreads; i++ ) {
		pthread_join(workerThreads[i], NULL);
	}
}
void WorkerPool::syncWorkerThreads() {
	pthread_mutex_lock(&jobsyncMutex);
	while ( jobsyncCounter > 0 ) {
		/*
		std::cout 	<< "\r$ waiting for " << jobsyncCounter 
					<< " jobs to finish..." << std::flush;
		*/
		pthread_cond_wait(&jobsyncConditionalVar, &jobsyncMutex);
	}
	pthread_mutex_unlock(&jobsyncMutex);
}
void WorkerPool::sendNoMoreJobsSignal() {
	pthread_mutex_lock(&jobsMutex);
	noMoreJobsFlag = true;
	pthread_cond_broadcast(&jobsConditionalVar); // wake all waiting threads
	pthread_mutex_unlock(&jobsMutex);
}
void WorkerPool::addNewJob(Job* job) {
	pthread_mutex_lock(&jobsyncMutex);
	jobsyncCounter++;
	pthread_mutex_unlock(&jobsyncMutex);

	pthread_mutex_lock(&jobsMutex);
	jobs.push_back(job);
	pthread_cond_signal(&jobsConditionalVar); // wake a waiting thread
	pthread_mutex_unlock(&jobsMutex);
}
void* WorkerPool::workerThreadEntryFunc(void* arg) {
	WorkerArg* warg = (WorkerArg*) arg; // pthread deals with void*s
	warg->wp->workerThreadEntry(warg->seed);
	return NULL;
}
void WorkerPool::workerThreadEntry(unsigned long seed) {
	// init a thread-local RNG
	Rng tlRNG(seed);
	// start the routine
	Job* myJob = NULL; bool lookForAJob = true;
	while ( lookForAJob ) {
		if ( myJob != NULL ) {
			bool ok = myJob->execute(tlRNG); // execute the job
			// todo: what to do with the return value?
			myJob->setFinished();
			if ( myJob->deleteWhenFinished ) delete myJob; // this is a bit dangerous...
			myJob = NULL; // reset my job
			pthread_mutex_lock(&jobsyncMutex);
			jobsyncCounter--;
			pthread_cond_signal(&jobsyncConditionalVar); 
			// syncWorkerThreads() could be listening... TODO: broadcast?
			pthread_mutex_unlock(&jobsyncMutex);
		}
		else {
			pthread_mutex_lock(&jobsMutex);
			if ( jobs.empty() ) {
				if ( noMoreJobsFlag ) lookForAJob = false;
				else {
					pthread_cond_wait(&jobsConditionalVar, &jobsMutex);
					if ( !jobs.empty() ) {
						myJob = jobs.front();
						jobs.pop_front();
					}
				}
			}
			else {
				myJob = jobs.front();
				jobs.pop_front();
			}
			pthread_mutex_unlock(&jobsMutex);
		}
	} // while ( lookForAJob )
}
