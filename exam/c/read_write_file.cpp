#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "ttas_sleep.lock.c"
#include <thread>
#include <iostream>
#include <fstream>
//using namespace std;


void * func_write_my_file(int arg, FILE * text_file){
    
    char of_name[10];
    sprintf(of_name,"test%d.txt",arg);
    //std::ofstream outp_file( "134.f" );
    FILE* outp_file = fopen(of_name, "w");
    char c;
    while(1){
        lock();
        while ( 1){
            //text_file>>c;
            //outp_file<<c;
            if (fscanf(text_file,"%c", &c)==EOF){
               fclose(outp_file);
               unlock();
               return 0;
            }
            if ( c == ' ' || c == '\n' || c == 0){
               fprintf(outp_file," "); 
               break;
            }
            fprintf(outp_file, "%c", c);
        }
        //outp_file<<std::endl;
        unlock();
    }
}

int main(int argc, char** argv)
{
        FILE *text_file;
        text_file = fopen("text.txt", "r");
        //std::ifstream text_file("text.txt");
	//long long time_start = __rdtsc();

	
	///* critical section
	std::thread thread_1(func_write_my_file, 1, text_file );
	std::thread thread_2(func_write_my_file, 2, text_file);
        std::thread thread_3(func_write_my_file, 3, text_file);
        std::thread thread_4(func_write_my_file, 4, text_file);
	thread_1.join();
        thread_2.join();
        thread_3.join();
        thread_4.join();
	//cout<<"time = "<<__rdtsc() - time_start << endl;
	return 0;
}

