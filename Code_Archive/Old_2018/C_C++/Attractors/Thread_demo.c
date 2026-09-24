#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<conio.h>
#include<sys/types.h>
#include<sys/fcntl.h>
#include<math.h>
#include<assert.h>
#include<pthread.h>


pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

void* f(void* param){

    pthread_mutex_lock(&mtx);
    printf("Thread no %d started \n", (int)param);
    pthread_mutex_unlock(&mtx);
    return NULL;
}

int main(){

    pthread_t t[10];
    int i = 0;
    for(i = 0; i<10; i++){
        pthread_create(&t[i], NULL, f, (void*)i);
    }
    for(i = 0; i<10; i++){
        pthread_join(t[i], NULL);
    }
    return 0;
}
