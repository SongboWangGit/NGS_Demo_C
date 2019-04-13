#include <unistd.h>
#include <cstdio>
#define FILE_NAME_SIZE 255

typedef struct {
    unsigned long size, resident, share, text, lib, data, dt;
}statm_t;

statm_t get_statm(pid_t pid);

void showMemStat(pid_t pid);
void displayStatm(pid_t pid, statm_t statm);
statm_t get_statm(pid_t pid)
{
    statm_t result = {0,0,0,0,0,0,0};

    char FILE_NAME[FILE_NAME_SIZE] = "/proc/self/statm";
    sprintf(FILE_NAME, "/proc/%d/statm", pid);

    FILE *fp = fopen(FILE_NAME, "r");
    fscanf(fp, "%lu %lu %lu %lu %lu %lu %lu", &result.size, &result.resident, &result.share, &result.text, &result.lib, &result.data, &result.dt);
    fclose(fp);
    return result;
}

void displayStatm(pid_t pid, statm_t statm)
{
    printf("Process %d Memory Use:\n", pid);
    printf("\t size \t resident \t share \t text \t lib \t data \t dt \n");
    printf("==========================================================\n");
    printf("\t %lu \t %lu \t\t %lu \t %lu \t %lu \t %lu \t %lu \n", statm.size, statm.resident, statm.share, statm.text, statm.lib, statm.data, statm.dt);
}

void showMemStat(pid_t pid)
{
    statm_t statm = get_statm(pid);
    displayStatm(pid, statm);
}
