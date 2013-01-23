#ifndef PRIMES_H_INCLUDED
#define PRIMES_H_INCLUDED

#define NUM_PRIMES 78498

//takes a number which will determine the maximum prime we want to look at
//and an int in which to store how many primes were filled
//returns an int array of those primes

int *loadPrimes(int uLimit, int *ret)
{
    int *primes = malloc(NUM_PRIMES*sizeof(int));
    FILE *pFile = fopen("primes.TXT", "r");
    char *nStr = malloc(10);
    int i;
    for(i = 0; i < NUM_PRIMES; i++)
    {
        fgets(nStr, 10, pFile);
        primes[i] = atoi(nStr);
        if(primes[i] > uLimit)
            break;
    }
    //primes = realloc(primes, i*sizeof(int));
    fclose(pFile);
    if(i < NUM_PRIMES)
        *ret = i;
    else *ret = NUM_PRIMES - 1;
    return primes;
}

#endif // PRIMES_H_INCLUDED
