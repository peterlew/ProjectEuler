#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include"useful.h"
#include"primes.h"

// PROBLEM #78

#define ARR_SIZE 100000

int *partitions;
int *pentagonals;

void fillPent()
{
    int i, j;
    j = 0;
    pentagonals[j++] = 0;
    for(i = 1; j < ARR_SIZE - 1; i++)
    {
        pentagonals[j++] = (3*i*i - i)/2;
        pentagonals[j++] = (3*i*i + i)/2;
    }
}

void clearParts()
{
    int i;
    for(i = 0; i < ARR_SIZE; i++)
        partitions[i] = 0;
}

int partsOf(int n)
{
    if(partitions[n])
        return partitions[n];
    if(n < 2)
        return 1;
    int curPent;
    int pentIndex = 1;
    int mult = 1;
    int sum = 0;
    while((curPent = pentagonals[pentIndex]) <= n)
    {
        sum += mult*partsOf(n - curPent);
        if(pentIndex % 2 == 0)
            mult *= -1;
        pentIndex++;
        sum = sum % 1000000;
    }
    partitions[n] = sum;
    return sum;
}

int main()
{

    partitions = calloc(ARR_SIZE, sizeof(int));
    pentagonals = calloc(ARR_SIZE, sizeof(int));

    clearParts();
    fillPent();

    int i;
    for(i = 1; i < ARR_SIZE; i++)
    {
        if(partsOf(i) == 0)
            break;
    }

    printf("Found parts div by 1000000 at %d\n", i);

    return 0;
}

/*

int **partitions;
FILE *outFile;

int pFind(int k, int n)
{
    if(k > n)
        return 0;
    if(k == n)
        return 1;
    fprintf(outFile, "partitions[ %d ] [ %d ] accessed\n", k, n);
    return partitions[k][n];
}

int main()
{

    outFile = fopen("output.txt", "w");
    int range = 30000;
    partitions = calloc(range, sizeof(int*));
    int i, j;
    for(i = 0; i < range; i++)
        if((partitions[i] = calloc(range, sizeof(int))) == NULL)
            {
                printf("Calloc failed!\n");
                printf("Fun fact: an int is %d bytes\n", sizeof(int));
                printf("You tried to allocate %lu bytes!\n", (unsigned long)range*(unsigned long)range*sizeof(int));
                exit(EXIT_FAILURE);
            }
    for(i = 0; i < range; i++)
        partitions[i][i] = 1;
    for(i = 2; i < range; i++)
    {
        for(j = i - 1; j > 0; j--)
        {
            partitions[j][i] = (pFind(j+1, i) + pFind(j, i-j)) % 1000000;
            if(j != 1 && j < i-j)
                free(partitions[j][i-j]);
        }
    }

    //int check = 5;
    //printf("Partitions of %d: %llu\n", check, partitions[1][check]);

    for(i = 2; i < range; i++)
    {
        if(partitions[1][i] == 0)
        {
            printf("Found! at %d\n", i);
            break;
        }
    }
    if(i == range)
    {
        printf("None found!\n");
        printf("But look at these numbers! (output.txt)\n");
        for(j = 2; j < range; j++)
            fprintf(outFile, "%d: %d\n", j, partitions[1][j], j);
    }

    return 0;
}

/* PROBLEM #77

int main()
{

    int partitionCount[101];
    int sopf[101];
    sopf[0] = partitionCount[0] = 0;
    sopf[1] = partitionCount[1] = 0;
    int i, j;
    int pCount;
    int div;
    for(i = 2; i < 100; i++)
    {
        pCount = 0;
        div = i;
        for(j = 2; div != 1; j++)
        {
            if(div % j == 0)
            {
                div = divideAway(div, j);
                pCount += j;
            }
        }
        sopf[i] = pCount;
    }

    int parts;
    for(i = 2; i < 100; i++)
    {
        parts = sopf[i];
        for(j = 1; j < i; j++)
            parts += sopf[j]*partitionCount[i-j];
        parts /= i;
        partitionCount[i] = parts;
    }

    for(i = 0; i < 100; i++)
        printf("pps in %d: %d\n", i, partitionCount[i]);

    return 0;
}

/* PROBLEM #72

int main()
{
    int check = 1000000;
    unsigned long long sum = 0;
    double *totients = malloc((check+1)*sizeof(double));
    totients[0] = 0;
    totients[1] = 0;
    int i;
    for(i = 2; i <= check; i++)
        totients[i] = i;
    int hMany;
    int *primes = loadPrimes(check, &hMany);
    int j, p;
    double pFact;
    for(i = 0; i <= hMany; i++)
    {
        p = primes[i];
        pFact = 1.0 - 1.0/(double)p;
        for(j = 1; p*j <= check; j++)
            totients[p*j] *= pFact;
    }
    for(i = 2; i <= check; i++)
    {
        sum += round(totients[i]);
        if(i%100 == 0)
            printf("Sum after %d: %llu\n", i, sum);
    }
    printf("Totient Sum For %d: %llu\n", check, sum);

    return 0;
}

/*

int divideAway(int n, int d)
{
    if(n % d == 0)
        return divideAway(n/d, d);
    else return n;
}

int totient(int n)
{
    double c, pProd;
    pProd = n;
    for(c = 2.0; c <= n; c++)
    {
        if(n % (int)c == 0)
        {
            pProd *= (1 - 1/c);
            n = divideAway(n, (int)c);
        }
    }
    return (int)pProd;
}

int mobius(int n)
{
    int i;
    if(n == 1)
        return 1;
    int nPrimeFacts = 0;
    for(i = 2; i <= n; i++)
    {
        if(n % i == 0)
        {
            n /= i;
            if(n % i == 0)
                return 0;
            nPrimeFacts++;
        }
    }
    return expt(-1, nPrimeFacts);
}

unsigned long long totientSum(int n)
{
    int i;
    unsigned long long sum = 0;
    for(i = 1; i <= n; i++)
        sum += mobius(i)*(n/i)*(1 + (n/i));
    return sum/2;
}

int main()
{
    int check = 1000000;
    printf("Totient Sum For %d: %u\n", check, totientSum(check) - 1);
    return 0;
}

/* PROBLEM #70

#define BIG_PRIME 17891

int myHash(char *s)
{
    int prod = 1;
    int i;
    for(i = 0; i < strlen(s); i++)
        prod = (prod * ((s[i] << 5) % BIG_PRIME)) % BIG_PRIME;
    return prod;
}

int myHashInt(int n)
{
    char *s = malloc(sizeof(int)*8+1);
    itoa(n, s, 10);
    return myHash(s);
}

int isPalindrome(int n1, int n2)
{
    return myHashInt(n1) == myHashInt(n2);
}

int main()
{

    int hMany;
    int *primes = loadPrimes(10000, &hMany);

    //int hMany = 25;

    int i, j;
    int p,q,prod, tot;
    double cRatio, bRatio;
    bRatio = 2.0;
    int bProd;
    for(i = hMany; i >= 0; i--)
    {
        p = primes[i];
        for(j = hMany; j >= 0; j--)
        {
            q = primes[j];
            prod = p*q;
            tot = (p-1)*(q-1);
            if(isPalindrome(prod, tot))
            {
                cRatio = (double)prod/(double)tot;
                printf("Number: %d\nTotient: %d\nRatio: %f\n", prod, tot, cRatio);
                if(cRatio < bRatio && prod < 10000000)
                {
                    bRatio = cRatio;
                    bProd = prod;
                }
            }
        }
    }

    printf("Best Number: %d\nRatio: %f\n", bProd, bRatio);

    return 0;
}

/* PROBLEM #66

int iterateChak(int a, int b, int k, int d)
{
    if(k == 1)
        return abs(a);
    else if(k == -1)
        return a*a + d*b*b;
    int m = 1;
    while((a + m*b)%k != 0)
        m++;
    int lastAbs, curAbs;
    curAbs = lastAbs = abs(m*m - d);
    while(1)
    {
        curAbs = abs(m*m - d);
        if(curAbs > lastAbs)
        {
            m--;
            while((a + m*b)%k != 0)
                m--;
            break;
        }
        else
        {
            lastAbs = curAbs;
            m++;
            while((a + m*b)%k != 0)
                m++;
        }
    }
    printf("New Triple: (%d, %d, %d) (D = %d) (m = %d)\n", (a*m+d*b)/abs(k),(a+b*m)/abs(k),(m*m-d)/k, d, m);
    iterateChak((a*m + d*b)/abs(k), (a + b*m)/abs(k), (m*m - d)/k, d);
}

int findMinX(int d)
{
    int a, k;
    a = ceil(sqrt(d));
    k = a*a - d;
    return iterateChak(a, 1, k, d);
}

int main()
{
    findMinX(5);
    return 0;
}

int maint()
{
    int bestD, bestX, result;
    bestD = 5;
    bestX = 9;
    int i;
    //findMinX(151);
    for(i = 182; i <= 1000; i++)
    {
        if(sqrt(i) == floor(sqrt(i)))
            continue;
        if((result = findMinX(i)) > bestX)
        {
            bestX = result;
            bestD = i;
        }
        printf("Checked %d : min X was %d\n", i, result);
    }
    printf("Best D = %d\n", bestD);
    printf("With X = %d\n", bestX);
    return 0;
}

/* PROBLEM #64

typedef struct valueSet{
    int front;
    int num;
    int root;
    int add;
} valueSet;

typedef struct listNode{
    valueSet *value;
    struct listNode *next;
} listNode;

typedef struct listHead{
    listNode *first;
    listNode *last;
} listHead;

valueSet *iterateSet(valueSet *set);

int listAdd(listHead *top, valueSet *set)
{
    listNode *newNode = malloc(sizeof(listNode));
    newNode->value = set;
    newNode->next = NULL;
    top->last->next = newNode;
    top->last = newNode;
    return 1;
}

int listFind(listHead *top, valueSet *set)
{
    listNode *cur = top->first;
    while(cur != NULL)
    {
        if(eqSet(cur->value, set))
        {
            int endCount = 0;
            while(cur != NULL)
            {
                endCount++;
                cur = cur->next;
            }
            return endCount;
        }
        else cur = cur->next;
    }
    return 0;
}

int eqSet(valueSet *a, valueSet *b)
{
    return (a->front == b->front && a->num == b->num &&
            a->root == b->root && a->add == b->add);
}

int printSet(valueSet *set)
{
    printf("%d + %d/sqrt(%d)%d\n", set->front, set->num, set->root, set->add);
    return 1;
}

int findPeriod(int root)
{
    int x = floor(sqrt(root));
    valueSet *firstSet = malloc(sizeof(valueSet));
    firstSet->front = x;
    firstSet->num = 1;
    firstSet->root = root;
    firstSet->add = x*-1;
    listNode *firstNode = malloc(sizeof(listNode));
    firstNode->value = firstSet;
    firstNode->next = NULL;
    listHead *seenList =  malloc(sizeof(listHead));
    seenList->first = firstNode;
    seenList->last = firstNode;
    valueSet *curSet = firstSet;
    while(1)
    {
        valueSet *nextSet = iterateSet(curSet);
        int result;
        if((result = listFind(seenList, nextSet)) != 0)
            return result;
        else
        {
            listAdd(seenList, nextSet);
            curSet = nextSet;
        }
    }
}

valueSet *iterateSet(valueSet *set)
{
    valueSet *newSet = malloc(sizeof(valueSet));
    int topRoot = set->root;
    int topAdd = set->add * -1;
    int den = set->root - set->add*set->add;
    den = den/set->num;
    newSet->front = floor((sqrt(topRoot) + topAdd)/den);
    newSet->num = den;
    newSet->root = set->root;
    newSet->add = topAdd - newSet->front*den;
    return newSet;
}

int main()
{
    int oddCount = 0;
    int i, period;
    for(i = 2; i < 10001; i++)
    {
        if(sqrt(i) != floor(sqrt(i)))
        {
            period = findPeriod(i);
            if(period % 2 == 1)
                oddCount++;
        }
    }
    printf("#Odd Periods: %d\n", oddCount);
    return 0;
}

/* PROBLEM #60

#define MAX_PRIME 2000
#define GROUP_SIZE 5

struct listNode{
    int value;
    struct listNode *nextNode;
};

struct listHead{
    int length;
    int prime;
    struct listNode *firstNode;
    struct listNode *lastNode;
};

typedef struct listNode listNode;
typedef struct listHead listHead;

int primes[MAX_PRIME];
listHead pairs[MAX_PRIME];

listNode *listISect(listNode *l1, listNode *l2);

int main()
{
    int i, j;
    j = 0;
    i = 3;
    while(j < MAX_PRIME)
    {
        if(isPrime(i))
            primes[j++] = i;
        i++;
    }
    printf("Primes Filled!\n(Highest Prime : %d)\n", primes[MAX_PRIME - 1]);
    for(i=0; i < MAX_PRIME; i++)
    {
        pairs[i].length = 0;
        pairs[i].prime = primes[i];
        pairs[i].firstNode = NULL;
        pairs[i].lastNode = NULL;
    }
    for(i=0; i < MAX_PRIME; i++)
    {
        for(j=0; j < MAX_PRIME; j++)
        {
            if(primeWith(primes[i], primes[j]))
            {
                pairs[i].length++;
                listNode *newNode = malloc(sizeof(listNode));
                newNode->value = primes[j];
                newNode->nextNode = NULL;
                if(pairs[i].firstNode == NULL)
                {
                    pairs[i].firstNode = newNode;
                    pairs[i].lastNode = newNode;
                }
                else
                {
                    pairs[i].lastNode->nextNode = newNode;
                    pairs[i].lastNode = newNode;
                }
            }
        }
    }
    printf("Filled! That didn't take long\n");
    for(i=0; i < MAX_PRIME; i++)
    {
        findNISect(i, pairs[i].firstNode, GROUP_SIZE - 1);
        //printf("%d\n", i);
    }
    return 0;
}

int findNISect(int index, listNode *list, int left)
{
    //printf("%d:%d  ", index, left);
    index++;
    while(index < MAX_PRIME)
    {
        listNode *iSect;
        iSect = listISect(list, pairs[index].firstNode);
        if(iSect == NULL)
            ;
        else
        {
            if(lenList(iSect) > GROUP_SIZE - 1)
            {
                if(left < 2)
                {
                    listNode *temp = iSect;
                    while(temp != NULL)
                    {
                        printf("%d\n", temp->value);
                        temp = temp->nextNode;
                    }
                    printf("DONE!\n");
                    exit(EXIT_SUCCESS);
                }
                else
                {
                    findNISect(index, iSect, left - 1);
                }
            }
        }
        index++;
    }
    //printf("Checked 'em all at ");
    return 0;
}

int lenList(listNode *list)
{
    int len = 0;
    while(list != NULL)
    {
        list = list->nextNode;
        len++;
    }
    return len;
}

listNode *listISect(listNode *l1, listNode *l2)
{
    listHead *newList = malloc(sizeof(listHead));
    newList->firstNode = NULL;
    newList->lastNode = NULL;
    while(l1 != NULL && l2 != NULL)
    {
        if(l1->value == l2->value)
        {
            //printf("common ground: %d\n", l1->value);
            listNode *newNode = malloc(sizeof(listNode));
            newNode->value = l1->value;
            newNode->nextNode = NULL;
            if(newList->firstNode == NULL)
                newList->firstNode = newNode;
            else newList->lastNode->nextNode = newNode;
            newList->lastNode = newNode;
            l1 = l1->nextNode;
            l2 = l2->nextNode;
        }
        else
        {
            if(l1->value > l2->value)
                l2 = l2->nextNode;
            else l1 = l1->nextNode;
        }
    }
    listNode *returnNode = newList->firstNode;
    free(newList);
    return returnNode;
}

int primeWith(int x, int y)
{
    if(x == y)
        return 1;
    if(isPrime(concatenate(x, y)) &&
       isPrime(concatenate(y, x)))
       return 1;
       else return 0;
}

int concatenate(int x, int y)
{
    return expt(10, magnitude(y))*x + y;
}

int magnitude(int x)
{
    int i;
    for(i = 1; expt(10, i) <= x; i++)
        ;
    return i;
}

/*

double marbFunc(double w, double b);

int main()
{
    double best = 0;
    double temp;
    double bestW, bestB;
    double w, b;
    for(w = 1; w < 50; w++)
        for(b = 0; b <= 50; b++)
        {
            if((temp = marbFunc(w, b)) > best)
            {
                best = temp;
                bestW = w;
                bestB = b;
            }
        }
    printf("Best:\tWhite: %f\n\tBlack: %f\n", bestW, bestB);
    printf("Probability of Survival: %f\n", best);
    return 0;
}

double marbFunc(double w, double b)
{
    double result = ((w/(w+b))+((50-w)/(100-w-b)))/2;
    return result;
}

/*
int fillPrimes()
{
    int i, j;
    j = 0;
    for(i = 2; j < MAX_PRIMES; i++)
        if(isPrime(i))
            primes[j++] = i;
    return 1;
}

int primeWithKnown(int x)
{
    if(primeWith(x, 3) &&
       primeWith(x, 7) &&
       primeWith(x, 109) &&
       primeWith(x, 673))
       return 1;
       else return 0;
}

int primeWithFive(int a, int b, int c, int d, int e)
{
    if(primeWith(primes[a], primes[b]) &&
       primeWith(primes[a], primes[c]) &&
       primeWith(primes[a], primes[d]) &&
       primeWith(primes[a], primes[e]) &&
       primeWith(primes[b], primes[c]) &&
       primeWith(primes[b], primes[d]) &&
       primeWith(primes[b], primes[e]) &&
       primeWith(primes[c], primes[d]) &&
       primeWith(primes[c], primes[e]) &&
       primeWith(primes[d], primes[e]))
       return 1;
       else return 0;
}

int primeWithFour(int a, int b, int c, int d)
{
    if(primeWith(primes[a], primes[b]) &&
       primeWith(primes[a], primes[c]) &&
       primeWith(primes[a], primes[d]) &&
       primeWith(primes[b], primes[c]) &&
       primeWith(primes[b], primes[d]) &&
       primeWith(primes[c], primes[d]))
       return 1;
       else return 0;
}

int primeWith(int x, int y)
{
    if(isPrime(concatenate(x, y)) &&
       isPrime(concatenate(y, x)))
       return 1;
       else return 0;
}

int concatenate(int x, int y)
{
    return expt(10, magnitude(y))*x + y;
}

int magnitude(int x)
{
    int i;
    for(i = 1; expt(10, i) <= x; i++)
        ;
    return i;
}

/* PROBLEM #61

int figNums[6][150];
int ordindex;

int main()
{
    fillFigs(figNums);
    int seen[6];
    int ordered[6];
    int i;
    for(i = 0; i < 6; i++)
    {
        seen[i] = 0;
        ordered[i] = 0;
    }
    ordindex = 0;
    search(1, seen, ordered);
    return 0;
}

int isCyclic(int x, int y)
{
    if(y == 0)
        return 0;
    if(x == 1)
        return 1;
    if((x % 100)/10 == y / 1000 && x % 10 == (y / 100) % 10)
        return 1;
    else return 0;
}

int search(int cur, int seen[], int ordered[])
{
    int i, j, allSeen;
    allSeen = 1;
    for(i = 0; i < 6; i++)
    {
        if(seen[i] == 0)
        {
            allSeen = 0;
            for(j = 0; j < 150; j++)
            {
                if(isCyclic(cur, figNums[i][j]))
                {
                    seen[i] = figNums[i][j];
                    ordered[ordindex++] = figNums[i][j];
                    search(figNums[i][j], seen, ordered);
                    seen[i] = 0;
                    ordered[--ordindex] = 0;
                }
            }
        }
    }
    if(allSeen && isCyclic(ordered[5], ordered[0]))
    {
        int sum;
        sum = 0;
        for(i = 0; i < 6; i++){
            printf("%d\t", ordered[i]);
            sum += ordered[i];
        }
        printf("Sum: %d\n", sum);
    }
    return 1;
}

int fd(int x)
{
    if(x >= 1000 && x < 10000)
        return x;
    else return 0;
}

int fillFigs(int a[])
{
    int i;
    for(i = 0; i < 150; i++)
    {
        figNums[0][i] = fd(i*(i+1)/2);
        figNums[1][i] = fd(i*i);
        figNums[2][i] = fd(i*(3*i-1)/2);
        figNums[3][i] = fd(i*(2*i-1));
        figNums[4][i] = fd(i*(5*i-3)/2);
        figNums[5][i] = fd(i*(3*i-2));
    }
    return 1;
}


/* PROBLEM #102

int coordinates[6000];

double getDegree(int x, int y);
double toDegree(double rads);

int main()
{
    inputArray(coordinates);
    int i, count;
    count = 0;
    i = 0;
    while(i < 5999)
    {
        if(inSpread(getDegree(coordinates[i++], coordinates[i++]),
                    getDegree(coordinates[i++], coordinates[i++]),
                    getDegree(coordinates[i++], coordinates[i++])))
            count++;
    }
    printf("Final tally: %d\n", count);
    return 0;
}

int inSpread(double a, double b, double c)
{
    double big, small, difference, mid;
    if(a > b)
    {
        big = a;
        small = b;
    }
    else
    {
        big = b;
        small = a;
    }
    if(big - small <= 180)
    {
        difference = big - small;
        mid = small + difference/2;
    }
    else
    {
        difference = small + 360 - big;
        mid = big + difference/2;
        if(mid >= 360)
            mid -= 360;
    }
    mid += 180;
    if(mid >= 360)
        mid -=360;
    if(c >= mid)
    {
        if(c - mid <= difference/2)
            return 1;
        else if(360 - c + mid <= difference/2)
            return 1;
        else return 0;
    }
    else
    {
        if(mid - c <= difference/2)
            return 1;
        else if(360 - mid + c <= difference/2)
            return 1;
        else return 0;
    }
}

double getDegree(int x, int y)
{
    if(y == 0)
    {
        if(x >= 0)
            return 0.0;
        else return 180.0;
    }
    if(x == 0)
    {
        if(y >= 0)
            return 90.0;
        else return 270.0;
    }
    double degs;
    degs = atan((double)y / (double)x);
    degs = toDegree(degs);
    if(degs < 0)
        degs += 360;
    if(x < 0)
        degs += 180;
    if(degs >= 360)
        degs -=360;
    return degs;
}

double toDegree(double rads)
{
    return rads*360/(2*M_PI);
}

int inputArray(int a[])
{
    FILE *fp;
    fp = fopen("triangles.txt", "r");
    int cur, i;
    char trash;
    for(i = 0; i < 6000; i++)
    {
        fscanf(fp, "%d", &cur);
        fscanf(fp, "%c", &trash);
        a[i] = cur;
    }
    fclose(fp);
    return 1;
}

/* PROBLEM #62

struct Entry
{
    unsigned long value;
    unsigned long original;
    int count;
};

typedef struct Entry Entry;

Entry* hashTable[10000];

unsigned long sortNum(unsigned long x);

int main()
{
    int i;
    for(i = 1; insertNum(expt_u(i, 3)) < 4; i++)
    ;
    return 0;
}

int numDigits(unsigned long x)
{
    int i;
    for(i = 0; x / expt_u(10, i); i++)
    ;
    return i;
}

int sortArray(int a[], int max)
{
    int i, j, temp;
    for(i = 0; i < max; i++)
    {
        for(j = i + 1; j < max; j++)
        {
            if(a[i] > a[j])
            {
                temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }
        }
    }
    return 1;
}

unsigned long sortNum(unsigned long x)
{
    int i, n;
    n = numDigits(x);
    int digits[n];
    for(i = 0; i < n; i++)
        digits[i] = (x / (expt_u(10, i))) % 10;
    sortArray(digits, n);
    unsigned long sum;
    sum = 0;
    for(i = 0; i < n; i++)
        sum += digits[i]*expt_u(10, i);
    return sum;
}

int getDigit(unsigned long x, int pos)
{
    return (x / expt_u(10, pos)) % 10;
}

int digHash(unsigned long x)
{
    unsigned long hashVal;
    hashVal = 1;
    int i;
    for(i = 0; expt_u(10, i) < x; i++)
    {
        hashVal *= (1779033703 + 2*getDigit(x, i));
    }
    hashVal = hashVal % 10000;
    return (int)hashVal;
}

int insertNum(unsigned long x)
{
    int hashVal;
    unsigned long sorted;
    sorted = sortNum(x);
    hashVal = digHash(sorted);
    int loopCount;
    loopCount = 1;
    while(hashTable[hashVal] != NULL && hashTable[hashVal]->value != sorted && loopCount < 10001)
    {
        if(hashVal == 9999)
            hashVal = 0;
        else hashVal++;
        if(++loopCount >= 10000)
            printf("Hash Table Full!");
    }
    if(hashTable[hashVal] == NULL)
    {
        Entry* newEntry = malloc(sizeof(Entry));
        newEntry->value = sorted;
        newEntry->original = x;
        newEntry->count = 1;
        hashTable[hashVal] = newEntry;
        printf("Inserted %u at Hash %d\n", x, hashVal);
        return 1;
    }
    else
    {
        hashTable[hashVal]->count++;
        printf("%u found to be permutation of %u\n", x, hashTable[hashVal]->original);
        return hashTable[hashVal]->count;
    }
}
/* PROBLEM #76

int ways[101][101];

int main()
{
    int i, j;
    for(i = 0; i < 101; i++)
    {
        ways[i][0] = 0;
        ways[i][1] = 1;
        for(j = 2; j < 101; j++)
        {
            ways[i][j] = -1;
        }
    }
    printf("Ways for 100: %d", findways(100, 99));
    return 0;
}

int findways(int goal, int cap)
{
    if(cap > goal)
        cap = goal;
    if(ways[goal][cap] != -1)
        return ways[goal][cap];
    if(cap == 1 || (goal == 1 && cap > 0))
    {
        ways[goal][cap] = 1;
        return 1;
    }
    int i;
    int count;
    count = 0;
    if(goal > cap)
        count += findways(goal - cap, cap);
    if(goal == cap)
        count += 1;
    count += findways(goal, cap - 1);
    ways[goal][cap] = count;
    return count;
}


/*

PROBLEM #81

#define MATRIX_SIZE 80

int matrix[MATRIX_SIZE * MATRIX_SIZE];
int pathLength[MATRIX_SIZE * MATRIX_SIZE];

int main()
{
    inputArray(matrix);
    int i;
    for(i = 0; i < MATRIX_SIZE * MATRIX_SIZE; i++)
        pathLength[i] = 0;
    pathLength[MATRIX_SIZE * MATRIX_SIZE - 1] = matrix[MATRIX_SIZE * MATRIX_SIZE - 1];
    printf("Max Path: %d\n", maxPath(0));
    return 0;
}

int maxPath(int index)
{
    if(pathLength[index])
        return pathLength[index];
    int mPath;
    if((index + 1) % MATRIX_SIZE == 0){
        mPath = matrix[index] + maxPath(index + MATRIX_SIZE);
        pathLength[index] = maxPath;
        return mPath;
    }
    if(index >= MATRIX_SIZE * MATRIX_SIZE - MATRIX_SIZE){
        mPath = matrix[index] + maxPath(index + 1);
        pathLength[index] = maxPath;
        return mPath;
    }
    int rSum, dSum;
    rSum = matrix[index] + maxPath(index + 1);
    dSum = matrix[index] + maxPath(index + MATRIX_SIZE);
    if (rSum < dSum)
        mPath = rSum;
        else mPath = dSum;
    pathLength[index] = mPath;
    return mPath;
}

int inputArray(int a[])
{
    FILE *fp;
    fp = fopen("matrix.txt", "r");
    int cur, i;
    char trash;
    for(i = 0; i < MATRIX_SIZE * MATRIX_SIZE; i++)
    {
        fscanf(fp, "%d", &cur);
        fscanf(fp, "%c", &trash);
        a[i] = cur;
    }
    fclose(fp);
    return 1;
}
/*

PROBLEM #99

int main()
{
    int pairs[2000];
    inputArray(pairs);
    double best, cur;
    best = 0.0;
    int i, bestLine;
    for(i = 0; i < 2000; i++){
        cur = log(pairs[i]) * pairs[++i];
        if(cur > best){
            best = cur;
            bestLine = (i + 1) / 2;
            printf("%d, %d : %f\n", pairs[i - 1], pairs[i], cur);
        }
    }
    printf("Best line: %d\n", bestLine);
    return 0;
}

int inputArray(int a[], int size)
{
    FILE *fp;
    fp = fopen("base_exp.txt", "r");
    int cur, i;
    char trash;
    for(i = 0; i < 2000; i++)
    {
        fscanf(fp, "%d", &cur);
        fscanf(fp, "%c", &trash);
        a[i] = cur;
    }
    fclose(fp);
    return 1;
}

/*

PROBLEM #58

int main()
{
    int cur, i, j;
    double count, primeCount;
    count = cur = 1;
    primeCount = 0;
    for(i = 2; primeCount / count >= .1 || i < 4; i+=2){
        for(j = 0; j < 4; j++){
            cur += i;
            if(isPrime(cur))
                primeCount++;
            count++;
        }
    }
    printf("%d\n", cur);
    return 0;
}

/*

PROBLEM #92

int seen[10000000];

int main()
{
    int i;
    for(i = 0; i < 10000000; i++)
        seen[i] = 0;
    seen[1] = 1;
    seen[89] = 89;
    int count;
    count = 0;
    for(i = 1; i < 10000000; i++)
        if(goesTo(i) == 89)
            count++;
    printf("%d\n", count);
    return 0;
}

int goesTo(int x)
{
    if(seen[x])
        return seen[x];
    int result;
    result = goesTo(sumSquareDigits(x));
    seen[x] = result;
    return result;
}

int sumSquareDigits(int x)
{
    int a[7];
    numToArray(x, a, 7);
    int i, sum;
    sum = 0;
    for(i = 0; i < 7; i++)
        if(a[i] != -1)
            sum += (a[i] * a[i]);
    return sum;
}

/*

PROBLEM #59

int cypherLength;

int main()
{
    int cypher[] = {79,59,12,2,79,35,8,28,20,2,3,68,8,9,68,45,0,12,9,67,68,4,7,5,23,27,1,21,79,85,78,79,85,71,38,10,71,27,12,2,79,6,2,8,13,9,1,13,9,8,68,19,7,1,71,56,11,21,11,68,6,3,22,2,14,0,30,79,1,31,6,23,19,10,0,73,79,44,2,79,19,6,28,68,16,6,16,15,79,35,8,11,72,71,14,10,3,79,12,2,79,19,6,28,68,32,0,0,73,79,86,71,39,1,71,24,5,20,79,13,9,79,16,15,10,68,5,10,3,14,1,10,14,1,3,71,24,13,19,7,68,32,0,0,73,79,87,71,39,1,71,12,22,2,14,16,2,11,68,2,25,1,21,22,16,15,6,10,0,79,16,15,10,22,2,79,13,20,65,68,41,0,16,15,6,10,0,79,1,31,6,23,19,28,68,19,7,5,19,79,12,2,79,0,14,11,10,64,27,68,10,14,15,2,65,68,83,79,40,14,9,1,71,6,16,20,10,8,1,79,19,6,28,68,14,1,68,15,6,9,75,79,5,9,11,68,19,7,13,20,79,8,14,9,1,71,8,13,17,10,23,71,3,13,0,7,16,71,27,11,71,10,18,2,29,29,8,1,1,73,79,81,71,59,12,2,79,8,14,8,12,19,79,23,15,6,10,2,28,68,19,7,22,8,26,3,15,79,16,15,10,68,3,14,22,12,1,1,20,28,72,71,14,10,3,79,16,15,10,68,3,14,22,12,1,1,20,28,68,4,14,10,71,1,1,17,10,22,71,10,28,19,6,10,0,26,13,20,7,68,14,27,74,71,89,68,32,0,0,71,28,1,9,27,68,45,0,12,9,79,16,15,10,68,37,14,20,19,6,23,19,79,83,71,27,11,71,27,1,11,3,68,2,25,1,21,22,11,9,10,68,6,13,11,18,27,68,19,7,1,71,3,13,0,7,16,71,28,11,71,27,12,6,27,68,2,25,1,21,22,11,9,10,68,10,6,3,15,27,68,5,10,8,14,10,18,2,79,6,2,12,5,18,28,1,71,0,2,71,7,13,20,79,16,2,28,16,14,2,11,9,22,74,71,87,68,45,0,12,9,79,12,14,2,23,2,3,2,71,24,5,20,79,10,8,27,68,19,7,1,71,3,13,0,7,16,92,79,12,2,79,19,6,28,68,8,1,8,30,79,5,71,24,13,19,1,1,20,28,68,19,0,68,19,7,1,71,3,13,0,7,16,73,79,93,71,59,12,2,79,11,9,10,68,16,7,11,71,6,23,71,27,12,2,79,16,21,26,1,71,3,13,0,7,16,75,79,19,15,0,68,0,6,18,2,28,68,11,6,3,15,27,68,19,0,68,2,25,1,21,22,11,9,10,72,71,24,5,20,79,3,8,6,10,0,79,16,8,79,7,8,2,1,71,6,10,19,0,68,19,7,1,71,24,11,21,3,0,73,79,85,87,79,38,18,27,68,6,3,16,15,0,17,0,7,68,19,7,1,71,24,11,21,3,0,71,24,5,20,79,9,6,11,1,71,27,12,21,0,17,0,7,68,15,6,9,75,79,16,15,10,68,16,0,22,11,11,68,3,6,0,9,72,16,71,29,1,4,0,3,9,6,30,2,79,12,14,2,68,16,7,1,9,79,12,2,79,7,6,2,1,73,79,85,86,79,33,17,10,10,71,6,10,71,7,13,20,79,11,16,1,68,11,14,10,3,79,5,9,11,68,6,2,11,9,8,68,15,6,23,71,0,19,9,79,20,2,0,20,11,10,72,71,7,1,71,24,5,20,79,10,8,27,68,6,12,7,2,31,16,2,11,74,71,94,86,71,45,17,19,79,16,8,79,5,11,3,68,16,7,11,71,13,1,11,6,1,17,10,0,71,7,13,10,79,5,9,11,68,6,12,7,2,31,16,2,11,68,15,6,9,75,79,12,2,79,3,6,25,1,71,27,12,2,79,22,14,8,12,19,79,16,8,79,6,2,12,11,10,10,68,4,7,13,11,11,22,2,1,68,8,9,68,32,0,0,73,79,85,84,79,48,15,10,29,71,14,22,2,79,22,2,13,11,21,1,69,71,59,12,14,28,68,14,28,68,9,0,16,71,14,68,23,7,29,20,6,7,6,3,68,5,6,22,19,7,68,21,10,23,18,3,16,14,1,3,71,9,22,8,2,68,15,26,9,6,1,68,23,14,23,20,6,11,9,79,11,21,79,20,11,14,10,75,79,16,15,6,23,71,29,1,5,6,22,19,7,68,4,0,9,2,28,68,1,29,11,10,79,35,8,11,74,86,91,68,52,0,68,19,7,1,71,56,11,21,11,68,5,10,7,6,2,1,71,7,17,10,14,10,71,14,10,3,79,8,14,25,1,3,79,12,2,29,1,71,0,10,71,10,5,21,27,12,71,14,9,8,1,3,71,26,23,73,79,44,2,79,19,6,28,68,1,26,8,11,79,11,1,79,17,9,9,5,14,3,13,9,8,68,11,0,18,2,79,5,9,11,68,1,14,13,19,7,2,18,3,10,2,28,23,73,79,37,9,11,68,16,10,68,15,14,18,2,79,23,2,10,10,71,7,13,20,79,3,11,0,22,30,67,68,19,7,1,71,8,8,8,29,29,71,0,2,71,27,12,2,79,11,9,3,29,71,60,11,9,79,11,1,79,16,15,10,68,33,14,16,15,10,22,73};
    cypherLength = sizeof(cypher) / sizeof(int);
    char x, y, z;
    for(x = 'a'; x <= 'z'; x++)
        for(y = 'a'; y <= 'z'; y++)
            for(z = 'a'; z <= 'z'; z++)
                {
                    int cypher[] = {79,59,12,2,79,35,8,28,20,2,3,68,8,9,68,45,0,12,9,67,68,4,7,5,23,27,1,21,79,85,78,79,85,71,38,10,71,27,12,2,79,6,2,8,13,9,1,13,9,8,68,19,7,1,71,56,11,21,11,68,6,3,22,2,14,0,30,79,1,31,6,23,19,10,0,73,79,44,2,79,19,6,28,68,16,6,16,15,79,35,8,11,72,71,14,10,3,79,12,2,79,19,6,28,68,32,0,0,73,79,86,71,39,1,71,24,5,20,79,13,9,79,16,15,10,68,5,10,3,14,1,10,14,1,3,71,24,13,19,7,68,32,0,0,73,79,87,71,39,1,71,12,22,2,14,16,2,11,68,2,25,1,21,22,16,15,6,10,0,79,16,15,10,22,2,79,13,20,65,68,41,0,16,15,6,10,0,79,1,31,6,23,19,28,68,19,7,5,19,79,12,2,79,0,14,11,10,64,27,68,10,14,15,2,65,68,83,79,40,14,9,1,71,6,16,20,10,8,1,79,19,6,28,68,14,1,68,15,6,9,75,79,5,9,11,68,19,7,13,20,79,8,14,9,1,71,8,13,17,10,23,71,3,13,0,7,16,71,27,11,71,10,18,2,29,29,8,1,1,73,79,81,71,59,12,2,79,8,14,8,12,19,79,23,15,6,10,2,28,68,19,7,22,8,26,3,15,79,16,15,10,68,3,14,22,12,1,1,20,28,72,71,14,10,3,79,16,15,10,68,3,14,22,12,1,1,20,28,68,4,14,10,71,1,1,17,10,22,71,10,28,19,6,10,0,26,13,20,7,68,14,27,74,71,89,68,32,0,0,71,28,1,9,27,68,45,0,12,9,79,16,15,10,68,37,14,20,19,6,23,19,79,83,71,27,11,71,27,1,11,3,68,2,25,1,21,22,11,9,10,68,6,13,11,18,27,68,19,7,1,71,3,13,0,7,16,71,28,11,71,27,12,6,27,68,2,25,1,21,22,11,9,10,68,10,6,3,15,27,68,5,10,8,14,10,18,2,79,6,2,12,5,18,28,1,71,0,2,71,7,13,20,79,16,2,28,16,14,2,11,9,22,74,71,87,68,45,0,12,9,79,12,14,2,23,2,3,2,71,24,5,20,79,10,8,27,68,19,7,1,71,3,13,0,7,16,92,79,12,2,79,19,6,28,68,8,1,8,30,79,5,71,24,13,19,1,1,20,28,68,19,0,68,19,7,1,71,3,13,0,7,16,73,79,93,71,59,12,2,79,11,9,10,68,16,7,11,71,6,23,71,27,12,2,79,16,21,26,1,71,3,13,0,7,16,75,79,19,15,0,68,0,6,18,2,28,68,11,6,3,15,27,68,19,0,68,2,25,1,21,22,11,9,10,72,71,24,5,20,79,3,8,6,10,0,79,16,8,79,7,8,2,1,71,6,10,19,0,68,19,7,1,71,24,11,21,3,0,73,79,85,87,79,38,18,27,68,6,3,16,15,0,17,0,7,68,19,7,1,71,24,11,21,3,0,71,24,5,20,79,9,6,11,1,71,27,12,21,0,17,0,7,68,15,6,9,75,79,16,15,10,68,16,0,22,11,11,68,3,6,0,9,72,16,71,29,1,4,0,3,9,6,30,2,79,12,14,2,68,16,7,1,9,79,12,2,79,7,6,2,1,73,79,85,86,79,33,17,10,10,71,6,10,71,7,13,20,79,11,16,1,68,11,14,10,3,79,5,9,11,68,6,2,11,9,8,68,15,6,23,71,0,19,9,79,20,2,0,20,11,10,72,71,7,1,71,24,5,20,79,10,8,27,68,6,12,7,2,31,16,2,11,74,71,94,86,71,45,17,19,79,16,8,79,5,11,3,68,16,7,11,71,13,1,11,6,1,17,10,0,71,7,13,10,79,5,9,11,68,6,12,7,2,31,16,2,11,68,15,6,9,75,79,12,2,79,3,6,25,1,71,27,12,2,79,22,14,8,12,19,79,16,8,79,6,2,12,11,10,10,68,4,7,13,11,11,22,2,1,68,8,9,68,32,0,0,73,79,85,84,79,48,15,10,29,71,14,22,2,79,22,2,13,11,21,1,69,71,59,12,14,28,68,14,28,68,9,0,16,71,14,68,23,7,29,20,6,7,6,3,68,5,6,22,19,7,68,21,10,23,18,3,16,14,1,3,71,9,22,8,2,68,15,26,9,6,1,68,23,14,23,20,6,11,9,79,11,21,79,20,11,14,10,75,79,16,15,6,23,71,29,1,5,6,22,19,7,68,4,0,9,2,28,68,1,29,11,10,79,35,8,11,74,86,91,68,52,0,68,19,7,1,71,56,11,21,11,68,5,10,7,6,2,1,71,7,17,10,14,10,71,14,10,3,79,8,14,25,1,3,79,12,2,29,1,71,0,10,71,10,5,21,27,12,71,14,9,8,1,3,71,26,23,73,79,44,2,79,19,6,28,68,1,26,8,11,79,11,1,79,17,9,9,5,14,3,13,9,8,68,11,0,18,2,79,5,9,11,68,1,14,13,19,7,2,18,3,10,2,28,23,73,79,37,9,11,68,16,10,68,15,14,18,2,79,23,2,10,10,71,7,13,20,79,3,11,0,22,30,67,68,19,7,1,71,8,8,8,29,29,71,0,2,71,27,12,2,79,11,9,3,29,71,60,11,9,79,11,1,79,16,15,10,68,33,14,16,15,10,22,73};
                    applyKey(cypher, x, y, z);
                    if(containsGospel(cypher)){
                        printChars(cypher);
                        printf("%c %c %c\n", x, y, z);
                        printf("%d %d %d\n", x, y, z);
                        printf("Sum: %d\n", sumArray(cypher));
                       }
                }
    return 0;
}

int sumArray(int cypher[])
{
    int sum, i;
    sum = 0;
    for(i = 0; i < cypherLength; i++)
        sum += cypher[i];
    return sum;
}

int containsThe(int cypher[])
{
    int i;
    for(i = 0; i < cypherLength; i++)
    {
        if(cypher[i] == 'e' &&
           cypher[i - 1] == 'h' &&
           cypher[i - 2] == 't')
            return 1;
    }
    return 0;
}

int containsGospel(int cypher[])
{
    int i;
    for(i = 0; i < cypherLength; i++)
    {
        if((cypher[i] == 'l') &&
           (cypher[i - 1] == 'e') &&
           (cypher[i - 2] == 'p') &&
           (cypher[i - 3] == 's') &&
           (cypher[i  - 4] == 'o') &&
           (cypher[i - 5] == 'G'))
           return 1;
    }
    return 0;
}

void applyKey(int a[], char x, char y, char z)
{
    int index, i;
    for(i = 0; i < cypherLength; i++){
        switch(index){
            case 0:
                a[i] = XOR(a[i], x);
                index++;
                break;
            case 1:
                a[i] = XOR(a[i], y);
                index++;
                break;
            case 2:
                a[i] = XOR(a[i], z);
                index = 0;
                break;
        }
    }
}

void printArray(int a[])
{
    int i;
    for(i = 7; i >= 0; i--)
        printf("%d", a[i]);
    printf("\n");
}

void printChars(int a[])
{
    int i;
    for(i = 0; i < cypherLength; i++)
        printf("%c", a[i]);
    printf("\n");
}

int XOR(int x, int y)
{
    int a[8];
    int b[8];
    toBinary(x, a);
    toBinary(y, b);
    int i, sum;
    sum = 0;
    for(i = 0; i < 8; i++)
        if((a[i] || b[i]) && !(a[i] && b[i]))
            sum += expt(2, i);
    return sum;
}

int toBinary(int x, int a[])
{
    int i;
    for(i = 0; i < 8; i++){
        if(x % 2 == 1)
            a[i] = 1;
            else a[i] = 0;
        x /= 2;
    }
}

/*

PROBLEM #44

#define MAX_PENT 10000

int main()
{
    int i, j;
    unsigned long pent, pent2;
    for(i = 1; i < MAX_PENT; i++){
        pent = i * ((3 * i) - 1) / 2;
        for(j = i + 1; j < MAX_PENT; j++){
            pent2 = j * ((3 * j) - 1) / 2;
            if(isPent(pent + pent2) && isPent(pent2 - pent))
                printf("%u\t%u\nD: %u\n", pent, pent2, pent2 - pent);
        }
    }
    return 0;
}

int isPent(unsigned long x)
{
    return isWholePosRoot(1.5, -0.5, -1.0 * (double)x);
}


/*

PROBLEM #46

#define MAX_RANGE 10000

int main()
{
    int i;
    for(i = 3; i < MAX_RANGE; i += 2){
        if(isPrime(i))
            continue;
        if(!isSquarePP(i)){
            printf("%d\n", i);
            return 0;
        }
    }
    printf("NOT FOUND!");
    return 0;
}

int isSquare(int x)
{
    double root;
    root = sqrt(x);
    if(root == round(root))
        return 1;
    return 0;
}

int isSquarePP(int x)
{
    int prime, diff;
    prime = 1;
    while((prime = nextPrime(prime)) < x){
        diff = x - prime;
        if(diff % 2 != 0)
            continue;
        if(isSquare(diff / 2))
            return prime;
    }
    return 0;
}

int nextPrime(int n)
{
    if (isPrime(n + 1))
        return n + 1;
    nextPrime(n + 1);
}
/*

PROBLEM #47

#define MAX_RANGE 300000

int main()
{
    printf("%d\n", hasNPF(2910, 4, 1));
    int i, streak;
    streak = 0;
    for(i = 1; i < MAX_RANGE; i++){
        if(hasNPF(i, 4, 1))
            streak++;
            else streak = 0;
        if(streak == 4){
            printf("%d\n", i - 3);
            return 0;
        }
    }
    printf("NOT FOUND!");
    return 0;
}

int hasNPF(int x, int n, int startPrime)
{
    if(n == 0)
        return 0;
    int prime;
    prime = startPrime;
    while((prime = nextPrime(prime)) <= sqrt(x))
        if (x % prime == 0)
            return (hasNPF(x / prime, n - 1, prime));
    if(n == 1)
        return 1;
        else return 0;
}

int nextPrime(int n)
{
    if (isPrime(n + 1))
        return n + 1;
    nextPrime(n + 1);
}

/*

PROBLEM #97

#define MAX_DIGITS 10

int main()
{
    int a[MAX_DIGITS];
    int i;
    for(i = 0; i < MAX_DIGITS; i++)
        a[i] = 0;
    a[0] = 2;
    powArray(a, 7830457, 2);
    mult10(a, 28433, 0, 0);
    printArrayNum10(a);
    printf("Plus one, of course!\n");
    return 1;
}

int add10(int a[], int b[])
{
    int i, sum;
    for(i = 0; i < MAX_DIGITS; i++)
    {
        sum = a[i] + b[i];
        b[i] = sum % 10;
        if (i == MAX_DIGITS - 1)
            return 1;
        if (sum > 9)
            b[i + 1] += sum / 10;
    }
    return 1;
}

int powArray(int a[], unsigned long pow, int x)
{
    pow--;
    for(pow; pow > 0; pow--)
        mult10(a, x, 0, 0);
    return 1;
}

int mult10(int a[], int x, int pos, int carry)
{
    int prod;
    prod = a[pos] * x;
    prod += carry;
    a[pos] = prod % 10;
    if (pos == MAX_DIGITS - 1)
        return 1;
    mult10(a, x, pos + 1, prod / 10);
}

int printArrayNum10(int x[])
{
    int i;
    for(i = MAX_DIGITS - 1; i >= 0; i--)
        printf("%d", x[i]);
    printf("\n");
    return 1;
}

/*

PROBLEM #33

int main()
{
    int i, j;
    double num1, num2, den1, den2, quo;
    for(i = 10; i < 100; i++){
        num1 = i / 10;
        num2 = i % 10;
        for(j = i + 1; j < 100; j++){
            if(i % 10 == 0 && j % 10 == 0)
                continue;
            den1 = j / 10;
            den2 = j % 10;
            quo = (double)i / (double)j;
            if((num2 == den2 && isNear(num1 / den1, quo)) ||
               (num1 == den2 && isNear(num2 / den1, quo)))
                printf("%d %d\n", i, j);
            if(den2 == 0)
                continue;
            if((num2 == den1 && isNear(num1 / den2, quo)) ||
               (num1 == den1 && isNear(num2 / den2, quo)))
                printf("%d %d\n", i, j);
        }
    }
    return 0;
}

int isNear(double a, double b)
{
    double diff;
    diff = a - b;
    if(diff < .0001 && diff > -.0001)
        return 1;
    return 0;
}

/*

PROBLEM #38

#define CAP 1000

int main()
{
    int solutions[4 * CAP];
    int a, b, i;
    double c;
    for(i = 0; i < 4 * CAP; i++)
        solutions[i] = 0;
    for(a = 1; a < CAP; a++){
        for(b = a; b < CAP; b++){
            c = sqrt(a*a + b*b);
            if (c == round(c))
                solutions[a + b + (int)c]++;
        }
    }
    int best, bestVal;
    best = 0;
    bestVal = 0;
    for(i = 0; i < CAP; i++)
        if(solutions[i] > bestVal){
            best = i;
            bestVal = solutions[i];
        }
    printf("Best value: %d\n", best);
    return 0;
}

/*

PROBLEM #37

#define MAX_DIGITS 10

int main()
{
    int a[MAX_DIGITS];
    int count, i, cur, sum;
    count = sum = 0;
    for(i = 11; count < 11; i++){
        if(!isPrime(i))
            continue;
        if(passesLeft(i, a) && passesRight(i, a)){
            count++;
            sum += i;
            printf("%d\n", i);
        }
    }
    printf("Sum: %d\n", sum);
    return 0;
}

int passesLeft(int x, int a[])
{
    int cur;
    cur = x;
    while(cur > 10){
        numToArray(cur, a);
        truncateLeft(a);
        cur = arrayToNum(a);
        if(!isPrime(cur))
            return 0;
    }
    return 1;
}

int passesRight(int x, int a[])
{
    int cur;
    cur = x;
    while(cur > 10){
        numToArray(cur, a);
        truncateRight(a);
        cur = arrayToNum(a);
        if(!isPrime(cur))
            return 0;
    }
    return 1;
}

int truncateLeft(int a[])
{
    int i;
    for(i = MAX_DIGITS - 1; a[i] == -1; i--)
    ;
    a[i] = -1;
    return 1;
}

int truncateRight(int a[])
{
    int i;
    for(i = 0; a[i] != -1; i++)
        a[i] = a[i + 1];
    return 1;
}

int isPrime(int x)
{
    if(x < 2)
        return 0;
    int i;
    for(i = 2; i <= sqrt(x); i++)
        if (x % i == 0)
            return 0;
    return 1;
}

int numToArray(int x, int a[])
{
    int i;
    for(i = 0; i < MAX_DIGITS; i++)
        a[i] = (x / (expt(10, i))) % 10;
    for(i = MAX_DIGITS - 1; a[i] == 0; i--)
        a[i] = -1;
    return 1;
}

int arrayToNum(int a[])
{
    int i, sum;
    sum = i = 0;
    while(i < MAX_DIGITS && a[i] != -1)
        sum += (a[i] * expt(10, i++));
    return sum;
}

int printArrayNum(int a[])
{
    int i;
    for(i = MAX_DIGITS - 1; i >= 0; i--)
        printf("%d", a[i]);
    printf("\n");
    return 1;
}

int expt(int x, int y)
{
    int prod;
    prod = 1;
    for(y; y > 0; y--)
        prod *= x;
    return prod;
}

/*

PROBLEM #31

unsigned long long waysToMake(int x, int bound);

unsigned long long waysMark[200][6];

int main()
{
    int i, j;
    for(i = 0; i < 200; i++)
        for(j = 0; j < 6; j++)
            waysMark[i][j] = 0;
    for(j = 0; j < 6; j++)
        waysMark[0][j] = 1;
    printf("%u\n", waysToMake(200, 5));
    return 0;
}

unsigned long long waysToMake(int x, int bound)
{
    if(waysMark[x][bound])
        return waysMark[x][bound];
    unsigned long long ways;
    ways = 0;
    if(x >= 1)
        ways += 1;
    if(x >= 2 && bound >= 0)
        ways += waysToMake(x - 2, 0);
    if(x >= 5 && bound >= 1)
        ways += waysToMake(x - 5, 1);
    if(x >= 10 && bound >= 2)
        ways += waysToMake(x - 10, 2);
    if(x >= 20 && bound >= 3)
        ways += waysToMake(x - 20, 3);
    if(x >= 50 && bound >= 4)
        ways += waysToMake(x - 50, 4);
    if(x >= 100 && bound >= 5)
        ways += waysToMake(x - 100, 5);
    if(x == 200)
        ways += 1;
    waysMark[x][bound] = ways;
    return ways;
}

/*

PROBLEM #52

#define MAX_DIGITS 10

int main()
{
    int a[MAX_DIGITS];
    int i, j, times, sorted;
    for(i = 1; i < 10000000; i++){
        numToArray(i, a);
        insertSort(a);
        sorted = arrayToNum(a);
        for(j = 2; j < 7; j++){
            times = j * i;
            numToArray(times, a);
            insertSort(a);
            if (sorted != arrayToNum(a))
                break;
        }
        if (j == 7)
            break;
    }
    printf("Number: %d\n", i);
    return 0;
}

int insertSort(int a[])
{
    int sorted[MAX_DIGITS];
    int max, i, j, jSaved;
    for(i = 0; i < MAX_DIGITS; i++){
        max = a[i];
        for(j = jSaved = i; j < MAX_DIGITS; j++){
            if(a[j] > max){
                max = a[j];
                jSaved = j;
            }
        }
        a[jSaved] = a[i];
        a[i] = max;
    }
    return 1;
}

int numToArray(int x, int a[])
{
    int i;
    for(i = 0; i < MAX_DIGITS; i++)
        a[i] = (x / (expt(10, i))) % 10;
    for(i = MAX_DIGITS - 1; a[i] == 0; i--)
        a[i] = -1;
    return 1;
}

int arrayToNum(int a[])
{
    int i, sum;
    sum = i = 0;
    while(i < MAX_DIGITS && a[i] != -1)
        sum += (a[i] * expt(10, i++));
    return sum;
}

int printArrayNum(int a[], int digits)
{
    int i;
    for(i = digits - 1; i >= 0; i--)
        printf("%d", a[i]);
    printf("\n");
    return 1;
}

int expt(int x, int y)
{
    int prod;
    prod = 1;
    for(y; y > 0; y--)
        prod *= x;
    return prod;
}

/*

PROBLEM #26

int main()
{
    int cur, best, bestIndex, i;
    best = 0;
    for(i = 2; i < 1000; i++){
        cur = cycleLength(i);
        if(cur > best){
            bestIndex = i;
            best = cur;
        }
    }
    printf("D Value: %d\n", bestIndex);
    printf("Cycle Length: %d\n", best);
    return 0;
}

int cycleLength(int x)
{
    int seen[1000];
    int i, j, length;
    for(i = 0; i < 1000; i++)
        seen[i] = 0;
    int divide;
    divide = 10;
    j = 1;
    while(!isIn(seen, divide)){
        seen[j++] = divide;
        if(x > divide)
            divide *= 10;
            else divide = 10 * (divide % x);
    }
    length = j - isIn(seen, divide);
    return length;
}

int isIn(int seen[], int divide)
{
    if(divide == 0)
        return 1;
    int i;
    for(i = 0; i < 1000; i++)
        if (seen[i] == divide)
            return i;
    return 0;
}

/*

PROBLEM #45

int main()
{
    int triangle, num;
    triangle = 40755;
    num = 285;
    num++;
    triangle += num;
    while(!(isWholePosRoot(1.5, -0.5, (double)(-1 * triangle)) &&
            isWholePosRoot(2.0, -1.0, (double)(-1 * triangle)))){
        num++;
        triangle += num;
    }
    printf("Number: %d\n", triangle);
    return 0;
}

int isWholePosRoot(double a, double b, double c)
{
    double root1, root2;
    root1 = (((-1 * b) + sqrt(pow(b, 2) - (4 * a * c))) / (2 * a));
    root2 = ((-1 * b) - sqrt(pow(b, 2) - (4 * a * c))) / (2 * a);
    if(root1 > 0 && root1 == round(root1) ||
       root2 > 0 && root2 == round(root2))
       return 1;
    return 0;
}

/*

PROBLEM #27

int main()
{
    int i, j, cur, best, bestA, bestB;
    best = 0;
    for(i = -999; i < 1000; i++){
        for(j = 0; j < 1000; j++){
            if(!isPrime(j))
                continue;
            cur = consecPrimes(i, j);
            if (cur > best){
                bestA = i;
                bestB = j;
                best = cur;
            }
        }
    }
    printf("A = %d\nB = %d\n", bestA, bestB);
    printf("Consecutive Primes: %d\n", best);
    printf("Product of coefficients: %d\n", (bestA *bestB));
    return 0;
}

int consecPrimes(int a, int b)
{
    int i;
    for(i = 0; isPrime((i * i) + (a * i) + b); i++)
    ;
    return i;
}

int isPrime(int x)
{
    if(x < 2)
        return 0;
    int i;
    for(i = 2; i <= sqrt(x); i++)
        if (x % i == 0)
            return 0;
    return 1;
}


/*

PROBLEM #42

int main()
{
    char *words[] = {"A","ABILITY","ABLE","ABOUT","ABOVE","ABSENCE","ABSOLUTELY","ACADEMIC","ACCEPT","ACCESS","ACCIDENT","ACCOMPANY","ACCORDING","ACCOUNT","ACHIEVE","ACHIEVEMENT","ACID","ACQUIRE","ACROSS","ACT","ACTION","ACTIVE","ACTIVITY","ACTUAL","ACTUALLY","ADD","ADDITION","ADDITIONAL","ADDRESS","ADMINISTRATION","ADMIT","ADOPT","ADULT","ADVANCE","ADVANTAGE","ADVICE","ADVISE","AFFAIR","AFFECT","AFFORD","AFRAID","AFTER","AFTERNOON","AFTERWARDS","AGAIN","AGAINST","AGE","AGENCY","AGENT","AGO","AGREE","AGREEMENT","AHEAD","AID","AIM","AIR","AIRCRAFT","ALL","ALLOW","ALMOST","ALONE","ALONG","ALREADY","ALRIGHT","ALSO","ALTERNATIVE","ALTHOUGH","ALWAYS","AMONG","AMONGST","AMOUNT","AN","ANALYSIS","ANCIENT","AND","ANIMAL","ANNOUNCE","ANNUAL","ANOTHER","ANSWER","ANY","ANYBODY","ANYONE","ANYTHING","ANYWAY","APART","APPARENT","APPARENTLY","APPEAL","APPEAR","APPEARANCE","APPLICATION","APPLY","APPOINT","APPOINTMENT","APPROACH","APPROPRIATE","APPROVE","AREA","ARGUE","ARGUMENT","ARISE","ARM","ARMY","AROUND","ARRANGE","ARRANGEMENT","ARRIVE","ART","ARTICLE","ARTIST","AS","ASK","ASPECT","ASSEMBLY","ASSESS","ASSESSMENT","ASSET","ASSOCIATE","ASSOCIATION","ASSUME","ASSUMPTION","AT","ATMOSPHERE","ATTACH","ATTACK","ATTEMPT","ATTEND","ATTENTION","ATTITUDE","ATTRACT","ATTRACTIVE","AUDIENCE","AUTHOR","AUTHORITY","AVAILABLE","AVERAGE","AVOID","AWARD","AWARE","AWAY","AYE","BABY","BACK","BACKGROUND","BAD","BAG","BALANCE","BALL","BAND","BANK","BAR","BASE","BASIC","BASIS","BATTLE","BE","BEAR","BEAT","BEAUTIFUL","BECAUSE","BECOME","BED","BEDROOM","BEFORE","BEGIN","BEGINNING","BEHAVIOUR","BEHIND","BELIEF","BELIEVE","BELONG","BELOW","BENEATH","BENEFIT","BESIDE","BEST","BETTER","BETWEEN","BEYOND","BIG","BILL","BIND","BIRD","BIRTH","BIT","BLACK","BLOCK","BLOOD","BLOODY","BLOW","BLUE","BOARD","BOAT","BODY","BONE","BOOK","BORDER","BOTH","BOTTLE","BOTTOM","BOX","BOY","BRAIN","BRANCH","BREAK","BREATH","BRIDGE","BRIEF","BRIGHT","BRING","BROAD","BROTHER","BUDGET","BUILD","BUILDING","BURN","BUS","BUSINESS","BUSY","BUT","BUY","BY","CABINET","CALL","CAMPAIGN","CAN","CANDIDATE","CAPABLE","CAPACITY","CAPITAL","CAR","CARD","CARE","CAREER","CAREFUL","CAREFULLY","CARRY","CASE","CASH","CAT","CATCH","CATEGORY","CAUSE","CELL","CENTRAL","CENTRE","CENTURY","CERTAIN","CERTAINLY","CHAIN","CHAIR","CHAIRMAN","CHALLENGE","CHANCE","CHANGE","CHANNEL","CHAPTER","CHARACTER","CHARACTERISTIC","CHARGE","CHEAP","CHECK","CHEMICAL","CHIEF","CHILD","CHOICE","CHOOSE","CHURCH","CIRCLE","CIRCUMSTANCE","CITIZEN","CITY","CIVIL","CLAIM","CLASS","CLEAN","CLEAR","CLEARLY","CLIENT","CLIMB","CLOSE","CLOSELY","CLOTHES","CLUB","COAL","CODE","COFFEE","COLD","COLLEAGUE","COLLECT","COLLECTION","COLLEGE","COLOUR","COMBINATION","COMBINE","COME","COMMENT","COMMERCIAL","COMMISSION","COMMIT","COMMITMENT","COMMITTEE","COMMON","COMMUNICATION","COMMUNITY","COMPANY","COMPARE","COMPARISON","COMPETITION","COMPLETE","COMPLETELY","COMPLEX","COMPONENT","COMPUTER","CONCENTRATE","CONCENTRATION","CONCEPT","CONCERN","CONCERNED","CONCLUDE","CONCLUSION","CONDITION","CONDUCT","CONFERENCE","CONFIDENCE","CONFIRM","CONFLICT","CONGRESS","CONNECT","CONNECTION","CONSEQUENCE","CONSERVATIVE","CONSIDER","CONSIDERABLE","CONSIDERATION","CONSIST","CONSTANT","CONSTRUCTION","CONSUMER","CONTACT","CONTAIN","CONTENT","CONTEXT","CONTINUE","CONTRACT","CONTRAST","CONTRIBUTE","CONTRIBUTION","CONTROL","CONVENTION","CONVERSATION","COPY","CORNER","CORPORATE","CORRECT","COS","COST","COULD","COUNCIL","COUNT","COUNTRY","COUNTY","COUPLE","COURSE","COURT","COVER","CREATE","CREATION","CREDIT","CRIME","CRIMINAL","CRISIS","CRITERION","CRITICAL","CRITICISM","CROSS","CROWD","CRY","CULTURAL","CULTURE","CUP","CURRENT","CURRENTLY","CURRICULUM","CUSTOMER","CUT","DAMAGE","DANGER","DANGEROUS","DARK","DATA","DATE","DAUGHTER","DAY","DEAD","DEAL","DEATH","DEBATE","DEBT","DECADE","DECIDE","DECISION","DECLARE","DEEP","DEFENCE","DEFENDANT","DEFINE","DEFINITION","DEGREE","DELIVER","DEMAND","DEMOCRATIC","DEMONSTRATE","DENY","DEPARTMENT","DEPEND","DEPUTY","DERIVE","DESCRIBE","DESCRIPTION","DESIGN","DESIRE","DESK","DESPITE","DESTROY","DETAIL","DETAILED","DETERMINE","DEVELOP","DEVELOPMENT","DEVICE","DIE","DIFFERENCE","DIFFERENT","DIFFICULT","DIFFICULTY","DINNER","DIRECT","DIRECTION","DIRECTLY","DIRECTOR","DISAPPEAR","DISCIPLINE","DISCOVER","DISCUSS","DISCUSSION","DISEASE","DISPLAY","DISTANCE","DISTINCTION","DISTRIBUTION","DISTRICT","DIVIDE","DIVISION","DO","DOCTOR","DOCUMENT","DOG","DOMESTIC","DOOR","DOUBLE","DOUBT","DOWN","DRAW","DRAWING","DREAM","DRESS","DRINK","DRIVE","DRIVER","DROP","DRUG","DRY","DUE","DURING","DUTY","EACH","EAR","EARLY","EARN","EARTH","EASILY","EAST","EASY","EAT","ECONOMIC","ECONOMY","EDGE","EDITOR","EDUCATION","EDUCATIONAL","EFFECT","EFFECTIVE","EFFECTIVELY","EFFORT","EGG","EITHER","ELDERLY","ELECTION","ELEMENT","ELSE","ELSEWHERE","EMERGE","EMPHASIS","EMPLOY","EMPLOYEE","EMPLOYER","EMPLOYMENT","EMPTY","ENABLE","ENCOURAGE","END","ENEMY","ENERGY","ENGINE","ENGINEERING","ENJOY","ENOUGH","ENSURE","ENTER","ENTERPRISE","ENTIRE","ENTIRELY","ENTITLE","ENTRY","ENVIRONMENT","ENVIRONMENTAL","EQUAL","EQUALLY","EQUIPMENT","ERROR","ESCAPE","ESPECIALLY","ESSENTIAL","ESTABLISH","ESTABLISHMENT","ESTATE","ESTIMATE","EVEN","EVENING","EVENT","EVENTUALLY","EVER","EVERY","EVERYBODY","EVERYONE","EVERYTHING","EVIDENCE","EXACTLY","EXAMINATION","EXAMINE","EXAMPLE","EXCELLENT","EXCEPT","EXCHANGE","EXECUTIVE","EXERCISE","EXHIBITION","EXIST","EXISTENCE","EXISTING","EXPECT","EXPECTATION","EXPENDITURE","EXPENSE","EXPENSIVE","EXPERIENCE","EXPERIMENT","EXPERT","EXPLAIN","EXPLANATION","EXPLORE","EXPRESS","EXPRESSION","EXTEND","EXTENT","EXTERNAL","EXTRA","EXTREMELY","EYE","FACE","FACILITY","FACT","FACTOR","FACTORY","FAIL","FAILURE","FAIR","FAIRLY","FAITH","FALL","FAMILIAR","FAMILY","FAMOUS","FAR","FARM","FARMER","FASHION","FAST","FATHER","FAVOUR","FEAR","FEATURE","FEE","FEEL","FEELING","FEMALE","FEW","FIELD","FIGHT","FIGURE","FILE","FILL","FILM","FINAL","FINALLY","FINANCE","FINANCIAL","FIND","FINDING","FINE","FINGER","FINISH","FIRE","FIRM","FIRST","FISH","FIT","FIX","FLAT","FLIGHT","FLOOR","FLOW","FLOWER","FLY","FOCUS","FOLLOW","FOLLOWING","FOOD","FOOT","FOOTBALL","FOR","FORCE","FOREIGN","FOREST","FORGET","FORM","FORMAL","FORMER","FORWARD","FOUNDATION","FREE","FREEDOM","FREQUENTLY","FRESH","FRIEND","FROM","FRONT","FRUIT","FUEL","FULL","FULLY","FUNCTION","FUND","FUNNY","FURTHER","FUTURE","GAIN","GAME","GARDEN","GAS","GATE","GATHER","GENERAL","GENERALLY","GENERATE","GENERATION","GENTLEMAN","GET","GIRL","GIVE","GLASS","GO","GOAL","GOD","GOLD","GOOD","GOVERNMENT","GRANT","GREAT","GREEN","GREY","GROUND","GROUP","GROW","GROWING","GROWTH","GUEST","GUIDE","GUN","HAIR","HALF","HALL","HAND","HANDLE","HANG","HAPPEN","HAPPY","HARD","HARDLY","HATE","HAVE","HE","HEAD","HEALTH","HEAR","HEART","HEAT","HEAVY","HELL","HELP","HENCE","HER","HERE","HERSELF","HIDE","HIGH","HIGHLY","HILL","HIM","HIMSELF","HIS","HISTORICAL","HISTORY","HIT","HOLD","HOLE","HOLIDAY","HOME","HOPE","HORSE","HOSPITAL","HOT","HOTEL","HOUR","HOUSE","HOUSEHOLD","HOUSING","HOW","HOWEVER","HUGE","HUMAN","HURT","HUSBAND","I","IDEA","IDENTIFY","IF","IGNORE","ILLUSTRATE","IMAGE","IMAGINE","IMMEDIATE","IMMEDIATELY","IMPACT","IMPLICATION","IMPLY","IMPORTANCE","IMPORTANT","IMPOSE","IMPOSSIBLE","IMPRESSION","IMPROVE","IMPROVEMENT","IN","INCIDENT","INCLUDE","INCLUDING","INCOME","INCREASE","INCREASED","INCREASINGLY","INDEED","INDEPENDENT","INDEX","INDICATE","INDIVIDUAL","INDUSTRIAL","INDUSTRY","INFLUENCE","INFORM","INFORMATION","INITIAL","INITIATIVE","INJURY","INSIDE","INSIST","INSTANCE","INSTEAD","INSTITUTE","INSTITUTION","INSTRUCTION","INSTRUMENT","INSURANCE","INTEND","INTENTION","INTEREST","INTERESTED","INTERESTING","INTERNAL","INTERNATIONAL","INTERPRETATION","INTERVIEW","INTO","INTRODUCE","INTRODUCTION","INVESTIGATE","INVESTIGATION","INVESTMENT","INVITE","INVOLVE","IRON","IS","ISLAND","ISSUE","IT","ITEM","ITS","ITSELF","JOB","JOIN","JOINT","JOURNEY","JUDGE","JUMP","JUST","JUSTICE","KEEP","KEY","KID","KILL","KIND","KING","KITCHEN","KNEE","KNOW","KNOWLEDGE","LABOUR","LACK","LADY","LAND","LANGUAGE","LARGE","LARGELY","LAST","LATE","LATER","LATTER","LAUGH","LAUNCH","LAW","LAWYER","LAY","LEAD","LEADER","LEADERSHIP","LEADING","LEAF","LEAGUE","LEAN","LEARN","LEAST","LEAVE","LEFT","LEG","LEGAL","LEGISLATION","LENGTH","LESS","LET","LETTER","LEVEL","LIABILITY","LIBERAL","LIBRARY","LIE","LIFE","LIFT","LIGHT","LIKE","LIKELY","LIMIT","LIMITED","LINE","LINK","LIP","LIST","LISTEN","LITERATURE","LITTLE","LIVE","LIVING","LOAN","LOCAL","LOCATION","LONG","LOOK","LORD","LOSE","LOSS","LOT","LOVE","LOVELY","LOW","LUNCH","MACHINE","MAGAZINE","MAIN","MAINLY","MAINTAIN","MAJOR","MAJORITY","MAKE","MALE","MAN","MANAGE","MANAGEMENT","MANAGER","MANNER","MANY","MAP","MARK","MARKET","MARRIAGE","MARRIED","MARRY","MASS","MASTER","MATCH","MATERIAL","MATTER","MAY","MAYBE","ME","MEAL","MEAN","MEANING","MEANS","MEANWHILE","MEASURE","MECHANISM","MEDIA","MEDICAL","MEET","MEETING","MEMBER","MEMBERSHIP","MEMORY","MENTAL","MENTION","MERELY","MESSAGE","METAL","METHOD","MIDDLE","MIGHT","MILE","MILITARY","MILK","MIND","MINE","MINISTER","MINISTRY","MINUTE","MISS","MISTAKE","MODEL","MODERN","MODULE","MOMENT","MONEY","MONTH","MORE","MORNING","MOST","MOTHER","MOTION","MOTOR","MOUNTAIN","MOUTH","MOVE","MOVEMENT","MUCH","MURDER","MUSEUM","MUSIC","MUST","MY","MYSELF","NAME","NARROW","NATION","NATIONAL","NATURAL","NATURE","NEAR","NEARLY","NECESSARILY","NECESSARY","NECK","NEED","NEGOTIATION","NEIGHBOUR","NEITHER","NETWORK","NEVER","NEVERTHELESS","NEW","NEWS","NEWSPAPER","NEXT","NICE","NIGHT","NO","NOBODY","NOD","NOISE","NONE","NOR","NORMAL","NORMALLY","NORTH","NORTHERN","NOSE","NOT","NOTE","NOTHING","NOTICE","NOTION","NOW","NUCLEAR","NUMBER","NURSE","OBJECT","OBJECTIVE","OBSERVATION","OBSERVE","OBTAIN","OBVIOUS","OBVIOUSLY","OCCASION","OCCUR","ODD","OF","OFF","OFFENCE","OFFER","OFFICE","OFFICER","OFFICIAL","OFTEN","OIL","OKAY","OLD","ON","ONCE","ONE","ONLY","ONTO","OPEN","OPERATE","OPERATION","OPINION","OPPORTUNITY","OPPOSITION","OPTION","OR","ORDER","ORDINARY","ORGANISATION","ORGANISE","ORGANIZATION","ORIGIN","ORIGINAL","OTHER","OTHERWISE","OUGHT","OUR","OURSELVES","OUT","OUTCOME","OUTPUT","OUTSIDE","OVER","OVERALL","OWN","OWNER","PACKAGE","PAGE","PAIN","PAINT","PAINTING","PAIR","PANEL","PAPER","PARENT","PARK","PARLIAMENT","PART","PARTICULAR","PARTICULARLY","PARTLY","PARTNER","PARTY","PASS","PASSAGE","PAST","PATH","PATIENT","PATTERN","PAY","PAYMENT","PEACE","PENSION","PEOPLE","PER","PERCENT","PERFECT","PERFORM","PERFORMANCE","PERHAPS","PERIOD","PERMANENT","PERSON","PERSONAL","PERSUADE","PHASE","PHONE","PHOTOGRAPH","PHYSICAL","PICK","PICTURE","PIECE","PLACE","PLAN","PLANNING","PLANT","PLASTIC","PLATE","PLAY","PLAYER","PLEASE","PLEASURE","PLENTY","PLUS","POCKET","POINT","POLICE","POLICY","POLITICAL","POLITICS","POOL","POOR","POPULAR","POPULATION","POSITION","POSITIVE","POSSIBILITY","POSSIBLE","POSSIBLY","POST","POTENTIAL","POUND","POWER","POWERFUL","PRACTICAL","PRACTICE","PREFER","PREPARE","PRESENCE","PRESENT","PRESIDENT","PRESS","PRESSURE","PRETTY","PREVENT","PREVIOUS","PREVIOUSLY","PRICE","PRIMARY","PRIME","PRINCIPLE","PRIORITY","PRISON","PRISONER","PRIVATE","PROBABLY","PROBLEM","PROCEDURE","PROCESS","PRODUCE","PRODUCT","PRODUCTION","PROFESSIONAL","PROFIT","PROGRAM","PROGRAMME","PROGRESS","PROJECT","PROMISE","PROMOTE","PROPER","PROPERLY","PROPERTY","PROPORTION","PROPOSE","PROPOSAL","PROSPECT","PROTECT","PROTECTION","PROVE","PROVIDE","PROVIDED","PROVISION","PUB","PUBLIC","PUBLICATION","PUBLISH","PULL","PUPIL","PURPOSE","PUSH","PUT","QUALITY","QUARTER","QUESTION","QUICK","QUICKLY","QUIET","QUITE","RACE","RADIO","RAILWAY","RAIN","RAISE","RANGE","RAPIDLY","RARE","RATE","RATHER","REACH","REACTION","READ","READER","READING","READY","REAL","REALISE","REALITY","REALIZE","REALLY","REASON","REASONABLE","RECALL","RECEIVE","RECENT","RECENTLY","RECOGNISE","RECOGNITION","RECOGNIZE","RECOMMEND","RECORD","RECOVER","RED","REDUCE","REDUCTION","REFER","REFERENCE","REFLECT","REFORM","REFUSE","REGARD","REGION","REGIONAL","REGULAR","REGULATION","REJECT","RELATE","RELATION","RELATIONSHIP","RELATIVE","RELATIVELY","RELEASE","RELEVANT","RELIEF","RELIGION","RELIGIOUS","RELY","REMAIN","REMEMBER","REMIND","REMOVE","REPEAT","REPLACE","REPLY","REPORT","REPRESENT","REPRESENTATION","REPRESENTATIVE","REQUEST","REQUIRE","REQUIREMENT","RESEARCH","RESOURCE","RESPECT","RESPOND","RESPONSE","RESPONSIBILITY","RESPONSIBLE","REST","RESTAURANT","RESULT","RETAIN","RETURN","REVEAL","REVENUE","REVIEW","REVOLUTION","RICH","RIDE","RIGHT","RING","RISE","RISK","RIVER","ROAD","ROCK","ROLE","ROLL","ROOF","ROOM","ROUND","ROUTE","ROW","ROYAL","RULE","RUN","RURAL","SAFE","SAFETY","SALE","SAME","SAMPLE","SATISFY","SAVE","SAY","SCALE","SCENE","SCHEME","SCHOOL","SCIENCE","SCIENTIFIC","SCIENTIST","SCORE","SCREEN","SEA","SEARCH","SEASON","SEAT","SECOND","SECONDARY","SECRETARY","SECTION","SECTOR","SECURE","SECURITY","SEE","SEEK","SEEM","SELECT","SELECTION","SELL","SEND","SENIOR","SENSE","SENTENCE","SEPARATE","SEQUENCE","SERIES","SERIOUS","SERIOUSLY","SERVANT","SERVE","SERVICE","SESSION","SET","SETTLE","SETTLEMENT","SEVERAL","SEVERE","SEX","SEXUAL","SHAKE","SHALL","SHAPE","SHARE","SHE","SHEET","SHIP","SHOE","SHOOT","SHOP","SHORT","SHOT","SHOULD","SHOULDER","SHOUT","SHOW","SHUT","SIDE","SIGHT","SIGN","SIGNAL","SIGNIFICANCE","SIGNIFICANT","SILENCE","SIMILAR","SIMPLE","SIMPLY","SINCE","SING","SINGLE","SIR","SISTER","SIT","SITE","SITUATION","SIZE","SKILL","SKIN","SKY","SLEEP","SLIGHTLY","SLIP","SLOW","SLOWLY","SMALL","SMILE","SO","SOCIAL","SOCIETY","SOFT","SOFTWARE","SOIL","SOLDIER","SOLICITOR","SOLUTION","SOME","SOMEBODY","SOMEONE","SOMETHING","SOMETIMES","SOMEWHAT","SOMEWHERE","SON","SONG","SOON","SORRY","SORT","SOUND","SOURCE","SOUTH","SOUTHERN","SPACE","SPEAK","SPEAKER","SPECIAL","SPECIES","SPECIFIC","SPEECH","SPEED","SPEND","SPIRIT","SPORT","SPOT","SPREAD","SPRING","STAFF","STAGE","STAND","STANDARD","STAR","START","STATE","STATEMENT","STATION","STATUS","STAY","STEAL","STEP","STICK","STILL","STOCK","STONE","STOP","STORE","STORY","STRAIGHT","STRANGE","STRATEGY","STREET","STRENGTH","STRIKE","STRONG","STRONGLY","STRUCTURE","STUDENT","STUDIO","STUDY","STUFF","STYLE","SUBJECT","SUBSTANTIAL","SUCCEED","SUCCESS","SUCCESSFUL","SUCH","SUDDENLY","SUFFER","SUFFICIENT","SUGGEST","SUGGESTION","SUITABLE","SUM","SUMMER","SUN","SUPPLY","SUPPORT","SUPPOSE","SURE","SURELY","SURFACE","SURPRISE","SURROUND","SURVEY","SURVIVE","SWITCH","SYSTEM","TABLE","TAKE","TALK","TALL","TAPE","TARGET","TASK","TAX","TEA","TEACH","TEACHER","TEACHING","TEAM","TEAR","TECHNICAL","TECHNIQUE","TECHNOLOGY","TELEPHONE","TELEVISION","TELL","TEMPERATURE","TEND","TERM","TERMS","TERRIBLE","TEST","TEXT","THAN","THANK","THANKS","THAT","THE","THEATRE","THEIR","THEM","THEME","THEMSELVES","THEN","THEORY","THERE","THEREFORE","THESE","THEY","THIN","THING","THINK","THIS","THOSE","THOUGH","THOUGHT","THREAT","THREATEN","THROUGH","THROUGHOUT","THROW","THUS","TICKET","TIME","TINY","TITLE","TO","TODAY","TOGETHER","TOMORROW","TONE","TONIGHT","TOO","TOOL","TOOTH","TOP","TOTAL","TOTALLY","TOUCH","TOUR","TOWARDS","TOWN","TRACK","TRADE","TRADITION","TRADITIONAL","TRAFFIC","TRAIN","TRAINING","TRANSFER","TRANSPORT","TRAVEL","TREAT","TREATMENT","TREATY","TREE","TREND","TRIAL","TRIP","TROOP","TROUBLE","TRUE","TRUST","TRUTH","TRY","TURN","TWICE","TYPE","TYPICAL","UNABLE","UNDER","UNDERSTAND","UNDERSTANDING","UNDERTAKE","UNEMPLOYMENT","UNFORTUNATELY","UNION","UNIT","UNITED","UNIVERSITY","UNLESS","UNLIKELY","UNTIL","UP","UPON","UPPER","URBAN","US","USE","USED","USEFUL","USER","USUAL","USUALLY","VALUE","VARIATION","VARIETY","VARIOUS","VARY","VAST","VEHICLE","VERSION","VERY","VIA","VICTIM","VICTORY","VIDEO","VIEW","VILLAGE","VIOLENCE","VISION","VISIT","VISITOR","VITAL","VOICE","VOLUME","VOTE","WAGE","WAIT","WALK","WALL","WANT","WAR","WARM","WARN","WASH","WATCH","WATER","WAVE","WAY","WE","WEAK","WEAPON","WEAR","WEATHER","WEEK","WEEKEND","WEIGHT","WELCOME","WELFARE","WELL","WEST","WESTERN","WHAT","WHATEVER","WHEN","WHERE","WHEREAS","WHETHER","WHICH","WHILE","WHILST","WHITE","WHO","WHOLE","WHOM","WHOSE","WHY","WIDE","WIDELY","WIFE","WILD","WILL","WIN","WIND","WINDOW","WINE","WING","WINNER","WINTER","WISH","WITH","WITHDRAW","WITHIN","WITHOUT","WOMAN","WONDER","WONDERFUL","WOOD","WORD","WORK","WORKER","WORKING","WORKS","WORLD","WORRY","WORTH","WOULD","WRITE","WRITER","WRITING","WRONG","YARD","YEAH","YEAR","YES","YESTERDAY","YET","YOU","YOUNG","YOUR","YOURSELF","YOUTH"};
    int triangles[40];
    fillTriangles(triangles);
    int count, numWords;
    count = 0;
    numWords = countWords(words);
    int i;
    for(i = 0; i < numWords; i++){
        if(isTriangle(words[i], triangles)){
            count++;
            printf("%s\n", words[i]);
        }
    }
    printf("Total: %d\n", count);
    return 0;
}

int isTriangle(char *word, int triangles[])
{
    int i, score;
    score = wordScore(word);
    for(i = 0; i < 40; i++)
        if(score == triangles[i])
            return 1;
    return 0;
}

int countWords(char *words[])
{
    int i;
    for(i = 0; strcmp(words[i], "YOUTH"); i++)
    ;
    return i + 1;
}

int wordScore(char *word)
{
    int i, sum;
    sum = 0;
    for(i = 0; i < strlen(word); i++)
        sum += (word[i] - 'A' + 1);
    return sum;
}

int fillTriangles(int triangles[])
{
    int i;
    triangles[0] = 0;
    for(i = 1; i < 40; i++)
        triangles[i] = triangles[i - 1] + i;
    return 1;
}

/*

PROBLEM #40

int main()
{
    int i, j, k, numDig, cur, count, prod;
    count = k = 0;
    prod = 1;
    for(i = 1; k < 7; i++){
        numDig = numDigits(i);
        for(j = numDig - 1; j >= 0 && k < 7; j--){
            cur = (i / (expt(10, j))) % 10;
            count++;
            if(count == expt(10, k)){
                printf("%d\n", cur);
                prod *= cur;
                k++;
            }
        }
    }
    printf("Product: %d\n", prod);
    return 0;
}

int numDigits(int x)
{
    int i;
    for(i = 0; x >= expt(10, i); i++)
    ;
    return i;
}

int expt(int x, int y)
{
    int prod;
    prod = 1;
    for(y; y > 0; y--)
        prod *= x;
    return prod;
}

/*

PROBLEM #35

int main()
{
    int cycleArray[6];
    int i, cycled, count;
    count = 0;
    for(i = 2; i < 1000000; i++){
        if (isPrime(i)){
            count++;
            numToArray(i, cycleArray);
            cycle(cycleArray);
            cycled = arrayToNum(cycleArray);
            if (cycled == i)
                    printf("%d\n", i);
            while(cycled != i){
                if (!isPrime(cycled)){
                    count--;
                    break;
                }
                cycle(cycleArray);
                cycled = arrayToNum(cycleArray);
                if (cycled == i)
                    printf("%d\n", i);
            }
        }
    }
    printf("Circular primes: %d\n", count);
    return 0;
}

int isPrime(int x)
{
    int i;
    for(i = 2; i <= sqrt(x); i++)
        if (x % i == 0)
            return 0;
    return 1;
}

int cycle(int x[])
{
    int i;
    for(i = 5; x[i] == -1 ; i--)
    ;
    int top;
    top = x[i];
    for(i; i >= 1; i--)
        x[i] = x[i - 1];
    x[0] = top;
    return 1;
}

int numToArray(int x, int a[])
{
    int i;
    for(i = 0; i < 6; i++)
        a[i] = (x / (expt(10, i))) % 10;
    for(i = 5; a[i] == 0; i--)
        a[i] = -1;
    return 1;
}

int arrayToNum(int a[])
{
    int i, sum;
    sum = i = 0;
    while(i < 6 && a[i] != -1)
        sum += (a[i] * expt(10, i++));
    return sum;
}

int printArrayNum(int a[], int digits)
{
    int i;
    for(i = digits - 1; i >= 0; i--)
        printf("%d", a[i]);
    printf("\n");
    return 1;
}

int expt(int x, int y)
{
    int prod;
    prod = 1;
    for(y; y > 0; y--)
        prod *= x;
    return prod;
}

/*

PROBLEM #23

#define MAX 28123

int main()
{
    int abundMark[MAX];
    int i, j, count;
    j = 0;
    count = 0;
    abundMark[0] = 0;
    for(i = 1; i < MAX; i++){
        if (sumFact(i) > i){
            abundMark[i] = 1;
            count++;
        }
            else abundMark[i] = 0;
    }
    int abunds[count];
    for(i = 0; i < MAX; i++)
        if (abundMark[i] == 1)
            abunds[j++] = i;
    int sumMark[MAX];
    for(i = 0; i < MAX; i++)
        sumMark[i] = 0;
    for(i = 0; i < count; i++)
        for(j = 0; j < count; j++)
            sumMark[abunds[i] + abunds[j]] = 1;
    unsigned long sum;
    sum = 0;
    for(i = 1; i < MAX; i++)
        if (sumMark[i] == 0)
            sum += i;
    printf("Sum: %u\n", sum);
    return 0;
}

int sumFact(int x)
{
    int i, sum;
    sum = 1;
    for(i = 2; i <= sqrt(x); i++){
        if (x % i == 0){
            sum += i;
            if (i != sqrt(x))
                sum += (x / i);
        }
    }
    return sum;
}

/*

PROBLEM #34

#define CAP 10000000

int facts[10];

int main()
{
    fillFact(facts);
    int a[7];
    int i;
    unsigned long sum;
    sum = 0;
    for(i = 1; i < CAP; i++){
        numToArray(i, a);
        if (i == sumFactDigs(a)){
            sum += i;
            printf("%d\n", i);
        }
    }
    printf("Sum: %u\n", sum);
    return 0;
}

int sumFactDigs(int a[])
{
    int i;
    for(i = 6; a[i] == 0; i--)
    ;
    int sum;
    sum = 0;
    for(i; i >= 0; i--)
        sum += facts[a[i]];
    return sum;
}

int numToArray(int x, int a[])
{
    int i;
    for(i = 0; i < 7; i++)
        a[i] = (x / ((int)(pow(10, (double)i)))) % 10;
    return 1;
}

int fillFact(int x[])
{
    int i;
    x[0] = 1;
    for(i = 1; i < 10; i++)
        x[i] = x[i - 1] * i;
    return 1;
}

/*

PROBLEM #36

#define CAP 1000000

int main()
{
    int bTen[6];
    int bTwo[20];
    unsigned long sum;
    sum = 0;
    int i;
    for(i = 1; i < CAP; i++){
        numToArray(i, bTen, 10, 6);
        numToArray(i, bTwo, 2, 20);
        if(isPalindrome(bTen, 6) && isPalindrome(bTwo, 20)){
            sum += i;
            printf("%d\n", i);
        }
    }
    printf("Sum: %u\n", sum);
    return 0;
}

int isPalindrome(int a[], int digits)
{
    int i, j;
    j = 0;
    for(i = digits - 1; a[i] == 0; i--)
    ;
    while(i >= j){
        if (a[i] != a[j])
            return 0;
        i--;
        j++;
    }
    return 1;
}

int numToArray(int x, int a[], int base, int digits)
{
    int i;
    for(i = 0; i < digits; i++)
        a[i] = (x / ((int)(pow((double)base, (double)i)))) % base;
    return 1;
}

int printArrayNum(int a[], int digits)
{
    int i;
    for(i = digits - 1; i >= 0; i--)
        printf("%d", a[i]);
    printf("\n");
    return 1;
}

/*

PROBLEM #29

int main()
{
    unsigned long total;
    total = 0;
    total += (81 * howMany(1));
    total += (4 * howMany(2));
    total += howMany(4);
    total += howMany(6);
    printf("Total: %u\n", total);
    return 0;
}

int howMany(int rank)
{
    int i, count, curRank;
    count = 0;
    for(i = 2; i <= 100 * rank; i++){
        if (i % 100 == 0)
            curRank = (i / 100);
            else curRank = (i / 100) + 1;
        if (divAnyOver(i, curRank, rank))
            count++;
    }
    return count;
}

int divAnyOver(int x, int init, int rank)
{
    for(init; init <= rank; init++)
        if (x % init == 0)
            return 1;
    return 0;
}

/*

PROBLEM  #24

unsigned long fact(int x);

int main()
{
    printf("%u\n", fact(9));
    return 0;
}

unsigned long fact(int x)
{
    unsigned long prod;
    prod = 1;
    for(x; x > 0; x--)
        prod *= x;
    return prod;
}

/*

PROBLEM #30

#define CAP 1000000

int main()
{
    int a[6];
    unsigned long sum;
    sum = 0;
    int i;
    for(i = 2; i < CAP; i++){
        numToArray(i, a);
        if (i == digitsToFifth(a)){
            sum += i;
            printf("%d\n", i);
        }
    }
    printf("Sum: %u\n", sum);
    return 0;
}

int digitsToFifth(int a[])
{
    int sum, i;
    sum = 0;
    for(i = 0; i < 6; i++)
        sum += pow(a[i], 5);
    return sum;
}

int numToArray(int x, int a[])
{
    int i;
    for(i = 0; i < 6; i++){
        a[0] = x % 10;
        a[1] = (x / 10) % 10;
        a[2] = (x / 100) % 10;
        a[3] = (x / 1000) % 10;
        a[4] = (x / 10000) % 10;
        a[5] = x / 100000;
    }
    return 1;
}

int printArrayNum(int x[])
{
    int i;
    for(i = 5; i >= 0; i--)
        printf("%d", x[i]);
    printf("\n");
    return 1;
}

/*

PROBLEM #28

#define SIZE 1001

int main()
{
    unsigned long sum;
    int distance, value, i;
    value = sum = 1;
    distance = 0;
    while(value < SIZE * SIZE){
        distance += 2;
        for(i = 0; i < 4; i++){
            value += distance;
            sum += value;
        }
    }
    printf("Sum of diagonals: %u\n", sum);
    return 0;
}

/*

PROBLEM #19

int main()
{
    int day, weekday, month, year;
    day = 1;
    weekday = 2;
    month = 1;
    year = 1900;
    int count;
    count = 0;
    while(!(day == 31 && month == 12 && year == 2000)){
        if (weekday == 1 && day == 1 && year > 1900)
            count++;
        if (weekday == 7)
            weekday = 1;
            else weekday++;
        if (day == monthQuota(month, year)){
            day = 1;
            if (month == 12){
                    month = 1;
                    year++;
            }
            else month++;
        }
        else day++;
    }
    printf("Sundays on 1st: %d\n", count);
    return 0;
}

int monthQuota(int month, int year)
{
    if (month == 2){
        if (year % 4 == 0 && year % 100 != 0)
            return 29;
            else return 28;
    }
    switch (month)
    {
        case 1:
        case 3:
        case 5:
        case 7:
        case 8:
        case 10:
        case 12:
            return 31;
        case 4:
        case 6:
        case 9:
        case 11:
            return 30;
    }
}

/*

PROBLEM #22

int numNames;

int main()
{
    char *names[] = {"MARY","PATRICIA","LINDA","BARBARA","ELIZABETH","JENNIFER","MARIA","SUSAN","MARGARET","DOROTHY","LISA","NANCY","KAREN","BETTY","HELEN","SANDRA","DONNA","CAROL","RUTH","SHARON","MICHELLE","LAURA","SARAH","KIMBERLY","DEBORAH","JESSICA","SHIRLEY","CYNTHIA","ANGELA","MELISSA","BRENDA","AMY","ANNA","REBECCA","VIRGINIA","KATHLEEN","PAMELA","MARTHA","DEBRA","AMANDA","STEPHANIE","CAROLYN","CHRISTINE","MARIE","JANET","CATHERINE","FRANCES","ANN","JOYCE","DIANE","ALICE","JULIE","HEATHER","TERESA","DORIS","GLORIA","EVELYN","JEAN","CHERYL","MILDRED","KATHERINE","JOAN","ASHLEY","JUDITH","ROSE","JANICE","KELLY","NICOLE","JUDY","CHRISTINA","KATHY","THERESA","BEVERLY","DENISE","TAMMY","IRENE","JANE","LORI","RACHEL","MARILYN","ANDREA","KATHRYN","LOUISE","SARA","ANNE","JACQUELINE","WANDA","BONNIE","JULIA","RUBY","LOIS","TINA","PHYLLIS","NORMA","PAULA","DIANA","ANNIE","LILLIAN","EMILY","ROBIN","PEGGY","CRYSTAL","GLADYS","RITA","DAWN","CONNIE","FLORENCE","TRACY","EDNA","TIFFANY","CARMEN","ROSA","CINDY","GRACE","WENDY","VICTORIA","EDITH","KIM","SHERRY","SYLVIA","JOSEPHINE","THELMA","SHANNON","SHEILA","ETHEL","ELLEN","ELAINE","MARJORIE","CARRIE","CHARLOTTE","MONICA","ESTHER","PAULINE","EMMA","JUANITA","ANITA","RHONDA","HAZEL","AMBER","EVA","DEBBIE","APRIL","LESLIE","CLARA","LUCILLE","JAMIE","JOANNE","ELEANOR","VALERIE","DANIELLE","MEGAN","ALICIA","SUZANNE","MICHELE","GAIL","BERTHA","DARLENE","VERONICA","JILL","ERIN","GERALDINE","LAUREN","CATHY","JOANN","LORRAINE","LYNN","SALLY","REGINA","ERICA","BEATRICE","DOLORES","BERNICE","AUDREY","YVONNE","ANNETTE","JUNE","SAMANTHA","MARION","DANA","STACY","ANA","RENEE","IDA","VIVIAN","ROBERTA","HOLLY","BRITTANY","MELANIE","LORETTA","YOLANDA","JEANETTE","LAURIE","KATIE","KRISTEN","VANESSA","ALMA","SUE","ELSIE","BETH","JEANNE","VICKI","CARLA","TARA","ROSEMARY","EILEEN","TERRI","GERTRUDE","LUCY","TONYA","ELLA","STACEY","WILMA","GINA","KRISTIN","JESSIE","NATALIE","AGNES","VERA","WILLIE","CHARLENE","BESSIE","DELORES","MELINDA","PEARL","ARLENE","MAUREEN","COLLEEN","ALLISON","TAMARA","JOY","GEORGIA","CONSTANCE","LILLIE","CLAUDIA","JACKIE","MARCIA","TANYA","NELLIE","MINNIE","MARLENE","HEIDI","GLENDA","LYDIA","VIOLA","COURTNEY","MARIAN","STELLA","CAROLINE","DORA","JO","VICKIE","MATTIE","TERRY","MAXINE","IRMA","MABEL","MARSHA","MYRTLE","LENA","CHRISTY","DEANNA","PATSY","HILDA","GWENDOLYN","JENNIE","NORA","MARGIE","NINA","CASSANDRA","LEAH","PENNY","KAY","PRISCILLA","NAOMI","CAROLE","BRANDY","OLGA","BILLIE","DIANNE","TRACEY","LEONA","JENNY","FELICIA","SONIA","MIRIAM","VELMA","BECKY","BOBBIE","VIOLET","KRISTINA","TONI","MISTY","MAE","SHELLY","DAISY","RAMONA","SHERRI","ERIKA","KATRINA","CLAIRE","LINDSEY","LINDSAY","GENEVA","GUADALUPE","BELINDA","MARGARITA","SHERYL","CORA","FAYE","ADA","NATASHA","SABRINA","ISABEL","MARGUERITE","HATTIE","HARRIET","MOLLY","CECILIA","KRISTI","BRANDI","BLANCHE","SANDY","ROSIE","JOANNA","IRIS","EUNICE","ANGIE","INEZ","LYNDA","MADELINE","AMELIA","ALBERTA","GENEVIEVE","MONIQUE","JODI","JANIE","MAGGIE","KAYLA","SONYA","JAN","LEE","KRISTINE","CANDACE","FANNIE","MARYANN","OPAL","ALISON","YVETTE","MELODY","LUZ","SUSIE","OLIVIA","FLORA","SHELLEY","KRISTY","MAMIE","LULA","LOLA","VERNA","BEULAH","ANTOINETTE","CANDICE","JUANA","JEANNETTE","PAM","KELLI","HANNAH","WHITNEY","BRIDGET","KARLA","CELIA","LATOYA","PATTY","SHELIA","GAYLE","DELLA","VICKY","LYNNE","SHERI","MARIANNE","KARA","JACQUELYN","ERMA","BLANCA","MYRA","LETICIA","PAT","KRISTA","ROXANNE","ANGELICA","JOHNNIE","ROBYN","FRANCIS","ADRIENNE","ROSALIE","ALEXANDRA","BROOKE","BETHANY","SADIE","BERNADETTE","TRACI","JODY","KENDRA","JASMINE","NICHOLE","RACHAEL","CHELSEA","MABLE","ERNESTINE","MURIEL","MARCELLA","ELENA","KRYSTAL","ANGELINA","NADINE","KARI","ESTELLE","DIANNA","PAULETTE","LORA","MONA","DOREEN","ROSEMARIE","ANGEL","DESIREE","ANTONIA","HOPE","GINGER","JANIS","BETSY","CHRISTIE","FREDA","MERCEDES","MEREDITH","LYNETTE","TERI","CRISTINA","EULA","LEIGH","MEGHAN","SOPHIA","ELOISE","ROCHELLE","GRETCHEN","CECELIA","RAQUEL","HENRIETTA","ALYSSA","JANA","KELLEY","GWEN","KERRY","JENNA","TRICIA","LAVERNE","OLIVE","ALEXIS","TASHA","SILVIA","ELVIRA","CASEY","DELIA","SOPHIE","KATE","PATTI","LORENA","KELLIE","SONJA","LILA","LANA","DARLA","MAY","MINDY","ESSIE","MANDY","LORENE","ELSA","JOSEFINA","JEANNIE","MIRANDA","DIXIE","LUCIA","MARTA","FAITH","LELA","JOHANNA","SHARI","CAMILLE","TAMI","SHAWNA","ELISA","EBONY","MELBA","ORA","NETTIE","TABITHA","OLLIE","JAIME","WINIFRED","KRISTIE","MARINA","ALISHA","AIMEE","RENA","MYRNA","MARLA","TAMMIE","LATASHA","BONITA","PATRICE","RONDA","SHERRIE","ADDIE","FRANCINE","DELORIS","STACIE","ADRIANA","CHERI","SHELBY","ABIGAIL","CELESTE","JEWEL","CARA","ADELE","REBEKAH","LUCINDA","DORTHY","CHRIS","EFFIE","TRINA","REBA","SHAWN","SALLIE","AURORA","LENORA","ETTA","LOTTIE","KERRI","TRISHA","NIKKI","ESTELLA","FRANCISCA","JOSIE","TRACIE","MARISSA","KARIN","BRITTNEY","JANELLE","LOURDES","LAUREL","HELENE","FERN","ELVA","CORINNE","KELSEY","INA","BETTIE","ELISABETH","AIDA","CAITLIN","INGRID","IVA","EUGENIA","CHRISTA","GOLDIE","CASSIE","MAUDE","JENIFER","THERESE","FRANKIE","DENA","LORNA","JANETTE","LATONYA","CANDY","MORGAN","CONSUELO","TAMIKA","ROSETTA","DEBORA","CHERIE","POLLY","DINA","JEWELL","FAY","JILLIAN","DOROTHEA","NELL","TRUDY","ESPERANZA","PATRICA","KIMBERLEY","SHANNA","HELENA","CAROLINA","CLEO","STEFANIE","ROSARIO","OLA","JANINE","MOLLIE","LUPE","ALISA","LOU","MARIBEL","SUSANNE","BETTE","SUSANA","ELISE","CECILE","ISABELLE","LESLEY","JOCELYN","PAIGE","JONI","RACHELLE","LEOLA","DAPHNE","ALTA","ESTER","PETRA","GRACIELA","IMOGENE","JOLENE","KEISHA","LACEY","GLENNA","GABRIELA","KERI","URSULA","LIZZIE","KIRSTEN","SHANA","ADELINE","MAYRA","JAYNE","JACLYN","GRACIE","SONDRA","CARMELA","MARISA","ROSALIND","CHARITY","TONIA","BEATRIZ","MARISOL","CLARICE","JEANINE","SHEENA","ANGELINE","FRIEDA","LILY","ROBBIE","SHAUNA","MILLIE","CLAUDETTE","CATHLEEN","ANGELIA","GABRIELLE","AUTUMN","KATHARINE","SUMMER","JODIE","STACI","LEA","CHRISTI","JIMMIE","JUSTINE","ELMA","LUELLA","MARGRET","DOMINIQUE","SOCORRO","RENE","MARTINA","MARGO","MAVIS","CALLIE","BOBBI","MARITZA","LUCILE","LEANNE","JEANNINE","DEANA","AILEEN","LORIE","LADONNA","WILLA","MANUELA","GALE","SELMA","DOLLY","SYBIL","ABBY","LARA","DALE","IVY","DEE","WINNIE","MARCY","LUISA","JERI","MAGDALENA","OFELIA","MEAGAN","AUDRA","MATILDA","LEILA","CORNELIA","BIANCA","SIMONE","BETTYE","RANDI","VIRGIE","LATISHA","BARBRA","GEORGINA","ELIZA","LEANN","BRIDGETTE","RHODA","HALEY","ADELA","NOLA","BERNADINE","FLOSSIE","ILA","GRETA","RUTHIE","NELDA","MINERVA","LILLY","TERRIE","LETHA","HILARY","ESTELA","VALARIE","BRIANNA","ROSALYN","EARLINE","CATALINA","AVA","MIA","CLARISSA","LIDIA","CORRINE","ALEXANDRIA","CONCEPCION","TIA","SHARRON","RAE","DONA","ERICKA","JAMI","ELNORA","CHANDRA","LENORE","NEVA","MARYLOU","MELISA","TABATHA","SERENA","AVIS","ALLIE","SOFIA","JEANIE","ODESSA","NANNIE","HARRIETT","LORAINE","PENELOPE","MILAGROS","EMILIA","BENITA","ALLYSON","ASHLEE","TANIA","TOMMIE","ESMERALDA","KARINA","EVE","PEARLIE","ZELMA","MALINDA","NOREEN","TAMEKA","SAUNDRA","HILLARY","AMIE","ALTHEA","ROSALINDA","JORDAN","LILIA","ALANA","GAY","CLARE","ALEJANDRA","ELINOR","MICHAEL","LORRIE","JERRI","DARCY","EARNESTINE","CARMELLA","TAYLOR","NOEMI","MARCIE","LIZA","ANNABELLE","LOUISA","EARLENE","MALLORY","CARLENE","NITA","SELENA","TANISHA","KATY","JULIANNE","JOHN","LAKISHA","EDWINA","MARICELA","MARGERY","KENYA","DOLLIE","ROXIE","ROSLYN","KATHRINE","NANETTE","CHARMAINE","LAVONNE","ILENE","KRIS","TAMMI","SUZETTE","CORINE","KAYE","JERRY","MERLE","CHRYSTAL","LINA","DEANNE","LILIAN","JULIANA","ALINE","LUANN","KASEY","MARYANNE","EVANGELINE","COLETTE","MELVA","LAWANDA","YESENIA","NADIA","MADGE","KATHIE","EDDIE","OPHELIA","VALERIA","NONA","MITZI","MARI","GEORGETTE","CLAUDINE","FRAN","ALISSA","ROSEANN","LAKEISHA","SUSANNA","REVA","DEIDRE","CHASITY","SHEREE","CARLY","JAMES","ELVIA","ALYCE","DEIRDRE","GENA","BRIANA","ARACELI","KATELYN","ROSANNE","WENDI","TESSA","BERTA","MARVA","IMELDA","MARIETTA","MARCI","LEONOR","ARLINE","SASHA","MADELYN","JANNA","JULIETTE","DEENA","AURELIA","JOSEFA","AUGUSTA","LILIANA","YOUNG","CHRISTIAN","LESSIE","AMALIA","SAVANNAH","ANASTASIA","VILMA","NATALIA","ROSELLA","LYNNETTE","CORINA","ALFREDA","LEANNA","CAREY","AMPARO","COLEEN","TAMRA","AISHA","WILDA","KARYN","CHERRY","QUEEN","MAURA","MAI","EVANGELINA","ROSANNA","HALLIE","ERNA","ENID","MARIANA","LACY","JULIET","JACKLYN","FREIDA","MADELEINE","MARA","HESTER","CATHRYN","LELIA","CASANDRA","BRIDGETT","ANGELITA","JANNIE","DIONNE","ANNMARIE","KATINA","BERYL","PHOEBE","MILLICENT","KATHERYN","DIANN","CARISSA","MARYELLEN","LIZ","LAURI","HELGA","GILDA","ADRIAN","RHEA","MARQUITA","HOLLIE","TISHA","TAMERA","ANGELIQUE","FRANCESCA","BRITNEY","KAITLIN","LOLITA","FLORINE","ROWENA","REYNA","TWILA","FANNY","JANELL","INES","CONCETTA","BERTIE","ALBA","BRIGITTE","ALYSON","VONDA","PANSY","ELBA","NOELLE","LETITIA","KITTY","DEANN","BRANDIE","LOUELLA","LETA","FELECIA","SHARLENE","LESA","BEVERLEY","ROBERT","ISABELLA","HERMINIA","TERRA","CELINA","TORI","OCTAVIA","JADE","DENICE","GERMAINE","SIERRA","MICHELL","CORTNEY","NELLY","DORETHA","SYDNEY","DEIDRA","MONIKA","LASHONDA","JUDI","CHELSEY","ANTIONETTE","MARGOT","BOBBY","ADELAIDE","NAN","LEEANN","ELISHA","DESSIE","LIBBY","KATHI","GAYLA","LATANYA","MINA","MELLISA","KIMBERLEE","JASMIN","RENAE","ZELDA","ELDA","MA","JUSTINA","GUSSIE","EMILIE","CAMILLA","ABBIE","ROCIO","KAITLYN","JESSE","EDYTHE","ASHLEIGH","SELINA","LAKESHA","GERI","ALLENE","PAMALA","MICHAELA","DAYNA","CARYN","ROSALIA","SUN","JACQULINE","REBECA","MARYBETH","KRYSTLE","IOLA","DOTTIE","BENNIE","BELLE","AUBREY","GRISELDA","ERNESTINA","ELIDA","ADRIANNE","DEMETRIA","DELMA","CHONG","JAQUELINE","DESTINY","ARLEEN","VIRGINA","RETHA","FATIMA","TILLIE","ELEANORE","CARI","TREVA","BIRDIE","WILHELMINA","ROSALEE","MAURINE","LATRICE","YONG","JENA","TARYN","ELIA","DEBBY","MAUDIE","JEANNA","DELILAH","CATRINA","SHONDA","HORTENCIA","THEODORA","TERESITA","ROBBIN","DANETTE","MARYJANE","FREDDIE","DELPHINE","BRIANNE","NILDA","DANNA","CINDI","BESS","IONA","HANNA","ARIEL","WINONA","VIDA","ROSITA","MARIANNA","WILLIAM","RACHEAL","GUILLERMINA","ELOISA","CELESTINE","CAREN","MALISSA","LONA","CHANTEL","SHELLIE","MARISELA","LEORA","AGATHA","SOLEDAD","MIGDALIA","IVETTE","CHRISTEN","ATHENA","JANEL","CHLOE","VEDA","PATTIE","TESSIE","TERA","MARILYNN","LUCRETIA","KARRIE","DINAH","DANIELA","ALECIA","ADELINA","VERNICE","SHIELA","PORTIA","MERRY","LASHAWN","DEVON","DARA","TAWANA","OMA","VERDA","CHRISTIN","ALENE","ZELLA","SANDI","RAFAELA","MAYA","KIRA","CANDIDA","ALVINA","SUZAN","SHAYLA","LYN","LETTIE","ALVA","SAMATHA","ORALIA","MATILDE","MADONNA","LARISSA","VESTA","RENITA","INDIA","DELOIS","SHANDA","PHILLIS","LORRI","ERLINDA","CRUZ","CATHRINE","BARB","ZOE","ISABELL","IONE","GISELA","CHARLIE","VALENCIA","ROXANNA","MAYME","KISHA","ELLIE","MELLISSA","DORRIS","DALIA","BELLA","ANNETTA","ZOILA","RETA","REINA","LAURETTA","KYLIE","CHRISTAL","PILAR","CHARLA","ELISSA","TIFFANI","TANA","PAULINA","LEOTA","BREANNA","JAYME","CARMEL","VERNELL","TOMASA","MANDI","DOMINGA","SANTA","MELODIE","LURA","ALEXA","TAMELA","RYAN","MIRNA","KERRIE","VENUS","NOEL","FELICITA","CRISTY","CARMELITA","BERNIECE","ANNEMARIE","TIARA","ROSEANNE","MISSY","CORI","ROXANA","PRICILLA","KRISTAL","JUNG","ELYSE","HAYDEE","ALETHA","BETTINA","MARGE","GILLIAN","FILOMENA","CHARLES","ZENAIDA","HARRIETTE","CARIDAD","VADA","UNA","ARETHA","PEARLINE","MARJORY","MARCELA","FLOR","EVETTE","ELOUISE","ALINA","TRINIDAD","DAVID","DAMARIS","CATHARINE","CARROLL","BELVA","NAKIA","MARLENA","LUANNE","LORINE","KARON","DORENE","DANITA","BRENNA","TATIANA","SAMMIE","LOUANN","LOREN","JULIANNA","ANDRIA","PHILOMENA","LUCILA","LEONORA","DOVIE","ROMONA","MIMI","JACQUELIN","GAYE","TONJA","MISTI","JOE","GENE","CHASTITY","STACIA","ROXANN","MICAELA","NIKITA","MEI","VELDA","MARLYS","JOHNNA","AURA","LAVERN","IVONNE","HAYLEY","NICKI","MAJORIE","HERLINDA","GEORGE","ALPHA","YADIRA","PERLA","GREGORIA","DANIEL","ANTONETTE","SHELLI","MOZELLE","MARIAH","JOELLE","CORDELIA","JOSETTE","CHIQUITA","TRISTA","LOUIS","LAQUITA","GEORGIANA","CANDI","SHANON","LONNIE","HILDEGARD","CECIL","VALENTINA","STEPHANY","MAGDA","KAROL","GERRY","GABRIELLA","TIANA","ROMA","RICHELLE","RAY","PRINCESS","OLETA","JACQUE","IDELLA","ALAINA","SUZANNA","JOVITA","BLAIR","TOSHA","RAVEN","NEREIDA","MARLYN","KYLA","JOSEPH","DELFINA","TENA","STEPHENIE","SABINA","NATHALIE","MARCELLE","GERTIE","DARLEEN","THEA","SHARONDA","SHANTEL","BELEN","VENESSA","ROSALINA","ONA","GENOVEVA","COREY","CLEMENTINE","ROSALBA","RENATE","RENATA","MI","IVORY","GEORGIANNA","FLOY","DORCAS","ARIANA","TYRA","THEDA","MARIAM","JULI","JESICA","DONNIE","VIKKI","VERLA","ROSELYN","MELVINA","JANNETTE","GINNY","DEBRAH","CORRIE","ASIA","VIOLETA","MYRTIS","LATRICIA","COLLETTE","CHARLEEN","ANISSA","VIVIANA","TWYLA","PRECIOUS","NEDRA","LATONIA","LAN","HELLEN","FABIOLA","ANNAMARIE","ADELL","SHARYN","CHANTAL","NIKI","MAUD","LIZETTE","LINDY","KIA","KESHA","JEANA","DANELLE","CHARLINE","CHANEL","CARROL","VALORIE","LIA","DORTHA","CRISTAL","SUNNY","LEONE","LEILANI","GERRI","DEBI","ANDRA","KESHIA","IMA","EULALIA","EASTER","DULCE","NATIVIDAD","LINNIE","KAMI","GEORGIE","CATINA","BROOK","ALDA","WINNIFRED","SHARLA","RUTHANN","MEAGHAN","MAGDALENE","LISSETTE","ADELAIDA","VENITA","TRENA","SHIRLENE","SHAMEKA","ELIZEBETH","DIAN","SHANTA","MICKEY","LATOSHA","CARLOTTA","WINDY","SOON","ROSINA","MARIANN","LEISA","JONNIE","DAWNA","CATHIE","BILLY","ASTRID","SIDNEY","LAUREEN","JANEEN","HOLLI","FAWN","VICKEY","TERESSA","SHANTE","RUBYE","MARCELINA","CHANDA","CARY","TERESE","SCARLETT","MARTY","MARNIE","LULU","LISETTE","JENIFFER","ELENOR","DORINDA","DONITA","CARMAN","BERNITA","ALTAGRACIA","ALETA","ADRIANNA","ZORAIDA","RONNIE","NICOLA","LYNDSEY","KENDALL","JANINA","CHRISSY","AMI","STARLA","PHYLIS","PHUONG","KYRA","CHARISSE","BLANCH","SANJUANITA","RONA","NANCI","MARILEE","MARANDA","CORY","BRIGETTE","SANJUANA","MARITA","KASSANDRA","JOYCELYN","IRA","FELIPA","CHELSIE","BONNY","MIREYA","LORENZA","KYONG","ILEANA","CANDELARIA","TONY","TOBY","SHERIE","OK","MARK","LUCIE","LEATRICE","LAKESHIA","GERDA","EDIE","BAMBI","MARYLIN","LAVON","HORTENSE","GARNET","EVIE","TRESSA","SHAYNA","LAVINA","KYUNG","JEANETTA","SHERRILL","SHARA","PHYLISS","MITTIE","ANABEL","ALESIA","THUY","TAWANDA","RICHARD","JOANIE","TIFFANIE","LASHANDA","KARISSA","ENRIQUETA","DARIA","DANIELLA","CORINNA","ALANNA","ABBEY","ROXANE","ROSEANNA","MAGNOLIA","LIDA","KYLE","JOELLEN","ERA","CORAL","CARLEEN","TRESA","PEGGIE","NOVELLA","NILA","MAYBELLE","JENELLE","CARINA","NOVA","MELINA","MARQUERITE","MARGARETTE","JOSEPHINA","EVONNE","DEVIN","CINTHIA","ALBINA","TOYA","TAWNYA","SHERITA","SANTOS","MYRIAM","LIZABETH","LISE","KEELY","JENNI","GISELLE","CHERYLE","ARDITH","ARDIS","ALESHA","ADRIANE","SHAINA","LINNEA","KAROLYN","HONG","FLORIDA","FELISHA","DORI","DARCI","ARTIE","ARMIDA","ZOLA","XIOMARA","VERGIE","SHAMIKA","NENA","NANNETTE","MAXIE","LOVIE","JEANE","JAIMIE","INGE","FARRAH","ELAINA","CAITLYN","STARR","FELICITAS","CHERLY","CARYL","YOLONDA","YASMIN","TEENA","PRUDENCE","PENNIE","NYDIA","MACKENZIE","ORPHA","MARVEL","LIZBETH","LAURETTE","JERRIE","HERMELINDA","CAROLEE","TIERRA","MIRIAN","META","MELONY","KORI","JENNETTE","JAMILA","ENA","ANH","YOSHIKO","SUSANNAH","SALINA","RHIANNON","JOLEEN","CRISTINE","ASHTON","ARACELY","TOMEKA","SHALONDA","MARTI","LACIE","KALA","JADA","ILSE","HAILEY","BRITTANI","ZONA","SYBLE","SHERRYL","RANDY","NIDIA","MARLO","KANDICE","KANDI","DEB","DEAN","AMERICA","ALYCIA","TOMMY","RONNA","NORENE","MERCY","JOSE","INGEBORG","GIOVANNA","GEMMA","CHRISTEL","AUDRY","ZORA","VITA","VAN","TRISH","STEPHAINE","SHIRLEE","SHANIKA","MELONIE","MAZIE","JAZMIN","INGA","HOA","HETTIE","GERALYN","FONDA","ESTRELLA","ADELLA","SU","SARITA","RINA","MILISSA","MARIBETH","GOLDA","EVON","ETHELYN","ENEDINA","CHERISE","CHANA","VELVA","TAWANNA","SADE","MIRTA","LI","KARIE","JACINTA","ELNA","DAVINA","CIERRA","ASHLIE","ALBERTHA","TANESHA","STEPHANI","NELLE","MINDI","LU","LORINDA","LARUE","FLORENE","DEMETRA","DEDRA","CIARA","CHANTELLE","ASHLY","SUZY","ROSALVA","NOELIA","LYDA","LEATHA","KRYSTYNA","KRISTAN","KARRI","DARLINE","DARCIE","CINDA","CHEYENNE","CHERRIE","AWILDA","ALMEDA","ROLANDA","LANETTE","JERILYN","GISELE","EVALYN","CYNDI","CLETA","CARIN","ZINA","ZENA","VELIA","TANIKA","PAUL","CHARISSA","THOMAS","TALIA","MARGARETE","LAVONDA","KAYLEE","KATHLENE","JONNA","IRENA","ILONA","IDALIA","CANDIS","CANDANCE","BRANDEE","ANITRA","ALIDA","SIGRID","NICOLETTE","MARYJO","LINETTE","HEDWIG","CHRISTIANA","CASSIDY","ALEXIA","TRESSIE","MODESTA","LUPITA","LITA","GLADIS","EVELIA","DAVIDA","CHERRI","CECILY","ASHELY","ANNABEL","AGUSTINA","WANITA","SHIRLY","ROSAURA","HULDA","EUN","BAILEY","YETTA","VERONA","THOMASINA","SIBYL","SHANNAN","MECHELLE","LUE","LEANDRA","LANI","KYLEE","KANDY","JOLYNN","FERNE","EBONI","CORENE","ALYSIA","ZULA","NADA","MOIRA","LYNDSAY","LORRETTA","JUAN","JAMMIE","HORTENSIA","GAYNELL","CAMERON","ADRIA","VINA","VICENTA","TANGELA","STEPHINE","NORINE","NELLA","LIANA","LESLEE","KIMBERELY","ILIANA","GLORY","FELICA","EMOGENE","ELFRIEDE","EDEN","EARTHA","CARMA","BEA","OCIE","MARRY","LENNIE","KIARA","JACALYN","CARLOTA","ARIELLE","YU","STAR","OTILIA","KIRSTIN","KACEY","JOHNETTA","JOEY","JOETTA","JERALDINE","JAUNITA","ELANA","DORTHEA","CAMI","AMADA","ADELIA","VERNITA","TAMAR","SIOBHAN","RENEA","RASHIDA","OUIDA","ODELL","NILSA","MERYL","KRISTYN","JULIETA","DANICA","BREANNE","AUREA","ANGLEA","SHERRON","ODETTE","MALIA","LORELEI","LIN","LEESA","KENNA","KATHLYN","FIONA","CHARLETTE","SUZIE","SHANTELL","SABRA","RACQUEL","MYONG","MIRA","MARTINE","LUCIENNE","LAVADA","JULIANN","JOHNIE","ELVERA","DELPHIA","CLAIR","CHRISTIANE","CHAROLETTE","CARRI","AUGUSTINE","ASHA","ANGELLA","PAOLA","NINFA","LEDA","LAI","EDA","SUNSHINE","STEFANI","SHANELL","PALMA","MACHELLE","LISSA","KECIA","KATHRYNE","KARLENE","JULISSA","JETTIE","JENNIFFER","HUI","CORRINA","CHRISTOPHER","CAROLANN","ALENA","TESS","ROSARIA","MYRTICE","MARYLEE","LIANE","KENYATTA","JUDIE","JANEY","IN","ELMIRA","ELDORA","DENNA","CRISTI","CATHI","ZAIDA","VONNIE","VIVA","VERNIE","ROSALINE","MARIELA","LUCIANA","LESLI","KARAN","FELICE","DENEEN","ADINA","WYNONA","TARSHA","SHERON","SHASTA","SHANITA","SHANI","SHANDRA","RANDA","PINKIE","PARIS","NELIDA","MARILOU","LYLA","LAURENE","LACI","JOI","JANENE","DOROTHA","DANIELE","DANI","CAROLYNN","CARLYN","BERENICE","AYESHA","ANNELIESE","ALETHEA","THERSA","TAMIKO","RUFINA","OLIVA","MOZELL","MARYLYN","MADISON","KRISTIAN","KATHYRN","KASANDRA","KANDACE","JANAE","GABRIEL","DOMENICA","DEBBRA","DANNIELLE","CHUN","BUFFY","BARBIE","ARCELIA","AJA","ZENOBIA","SHAREN","SHAREE","PATRICK","PAGE","MY","LAVINIA","KUM","KACIE","JACKELINE","HUONG","FELISA","EMELIA","ELEANORA","CYTHIA","CRISTIN","CLYDE","CLARIBEL","CARON","ANASTACIA","ZULMA","ZANDRA","YOKO","TENISHA","SUSANN","SHERILYN","SHAY","SHAWANDA","SABINE","ROMANA","MATHILDA","LINSEY","KEIKO","JOANA","ISELA","GRETTA","GEORGETTA","EUGENIE","DUSTY","DESIRAE","DELORA","CORAZON","ANTONINA","ANIKA","WILLENE","TRACEE","TAMATHA","REGAN","NICHELLE","MICKIE","MAEGAN","LUANA","LANITA","KELSIE","EDELMIRA","BREE","AFTON","TEODORA","TAMIE","SHENA","MEG","LINH","KELI","KACI","DANYELLE","BRITT","ARLETTE","ALBERTINE","ADELLE","TIFFINY","STORMY","SIMONA","NUMBERS","NICOLASA","NICHOL","NIA","NAKISHA","MEE","MAIRA","LOREEN","KIZZY","JOHNNY","JAY","FALLON","CHRISTENE","BOBBYE","ANTHONY","YING","VINCENZA","TANJA","RUBIE","RONI","QUEENIE","MARGARETT","KIMBERLI","IRMGARD","IDELL","HILMA","EVELINA","ESTA","EMILEE","DENNISE","DANIA","CARL","CARIE","ANTONIO","WAI","SANG","RISA","RIKKI","PARTICIA","MUI","MASAKO","MARIO","LUVENIA","LOREE","LONI","LIEN","KEVIN","GIGI","FLORENCIA","DORIAN","DENITA","DALLAS","CHI","BILLYE","ALEXANDER","TOMIKA","SHARITA","RANA","NIKOLE","NEOMA","MARGARITE","MADALYN","LUCINA","LAILA","KALI","JENETTE","GABRIELE","EVELYNE","ELENORA","CLEMENTINA","ALEJANDRINA","ZULEMA","VIOLETTE","VANNESSA","THRESA","RETTA","PIA","PATIENCE","NOELLA","NICKIE","JONELL","DELTA","CHUNG","CHAYA","CAMELIA","BETHEL","ANYA","ANDREW","THANH","SUZANN","SPRING","SHU","MILA","LILLA","LAVERNA","KEESHA","KATTIE","GIA","GEORGENE","EVELINE","ESTELL","ELIZBETH","VIVIENNE","VALLIE","TRUDIE","STEPHANE","MICHEL","MAGALY","MADIE","KENYETTA","KARREN","JANETTA","HERMINE","HARMONY","DRUCILLA","DEBBI","CELESTINA","CANDIE","BRITNI","BECKIE","AMINA","ZITA","YUN","YOLANDE","VIVIEN","VERNETTA","TRUDI","SOMMER","PEARLE","PATRINA","OSSIE","NICOLLE","LOYCE","LETTY","LARISA","KATHARINA","JOSELYN","JONELLE","JENELL","IESHA","HEIDE","FLORINDA","FLORENTINA","FLO","ELODIA","DORINE","BRUNILDA","BRIGID","ASHLI","ARDELLA","TWANA","THU","TARAH","SUNG","SHEA","SHAVON","SHANE","SERINA","RAYNA","RAMONITA","NGA","MARGURITE","LUCRECIA","KOURTNEY","KATI","JESUS","JESENIA","DIAMOND","CRISTA","AYANA","ALICA","ALIA","VINNIE","SUELLEN","ROMELIA","RACHELL","PIPER","OLYMPIA","MICHIKO","KATHALEEN","JOLIE","JESSI","JANESSA","HANA","HA","ELEASE","CARLETTA","BRITANY","SHONA","SALOME","ROSAMOND","REGENA","RAINA","NGOC","NELIA","LOUVENIA","LESIA","LATRINA","LATICIA","LARHONDA","JINA","JACKI","HOLLIS","HOLLEY","EMMY","DEEANN","CORETTA","ARNETTA","VELVET","THALIA","SHANICE","NETA","MIKKI","MICKI","LONNA","LEANA","LASHUNDA","KILEY","JOYE","JACQULYN","IGNACIA","HYUN","HIROKO","HENRY","HENRIETTE","ELAYNE","DELINDA","DARNELL","DAHLIA","COREEN","CONSUELA","CONCHITA","CELINE","BABETTE","AYANNA","ANETTE","ALBERTINA","SKYE","SHAWNEE","SHANEKA","QUIANA","PAMELIA","MIN","MERRI","MERLENE","MARGIT","KIESHA","KIERA","KAYLENE","JODEE","JENISE","ERLENE","EMMIE","ELSE","DARYL","DALILA","DAISEY","CODY","CASIE","BELIA","BABARA","VERSIE","VANESA","SHELBA","SHAWNDA","SAM","NORMAN","NIKIA","NAOMA","MARNA","MARGERET","MADALINE","LAWANA","KINDRA","JUTTA","JAZMINE","JANETT","HANNELORE","GLENDORA","GERTRUD","GARNETT","FREEDA","FREDERICA","FLORANCE","FLAVIA","DENNIS","CARLINE","BEVERLEE","ANJANETTE","VALDA","TRINITY","TAMALA","STEVIE","SHONNA","SHA","SARINA","ONEIDA","MICAH","MERILYN","MARLEEN","LURLINE","LENNA","KATHERIN","JIN","JENI","HAE","GRACIA","GLADY","FARAH","ERIC","ENOLA","EMA","DOMINQUE","DEVONA","DELANA","CECILA","CAPRICE","ALYSHA","ALI","ALETHIA","VENA","THERESIA","TAWNY","SONG","SHAKIRA","SAMARA","SACHIKO","RACHELE","PAMELLA","NICKY","MARNI","MARIEL","MAREN","MALISA","LIGIA","LERA","LATORIA","LARAE","KIMBER","KATHERN","KAREY","JENNEFER","JANETH","HALINA","FREDIA","DELISA","DEBROAH","CIERA","CHIN","ANGELIKA","ANDREE","ALTHA","YEN","VIVAN","TERRESA","TANNA","SUK","SUDIE","SOO","SIGNE","SALENA","RONNI","REBBECCA","MYRTIE","MCKENZIE","MALIKA","MAIDA","LOAN","LEONARDA","KAYLEIGH","FRANCE","ETHYL","ELLYN","DAYLE","CAMMIE","BRITTNI","BIRGIT","AVELINA","ASUNCION","ARIANNA","AKIKO","VENICE","TYESHA","TONIE","TIESHA","TAKISHA","STEFFANIE","SINDY","SANTANA","MEGHANN","MANDA","MACIE","LADY","KELLYE","KELLEE","JOSLYN","JASON","INGER","INDIRA","GLINDA","GLENNIS","FERNANDA","FAUSTINA","ENEIDA","ELICIA","DOT","DIGNA","DELL","ARLETTA","ANDRE","WILLIA","TAMMARA","TABETHA","SHERRELL","SARI","REFUGIO","REBBECA","PAULETTA","NIEVES","NATOSHA","NAKITA","MAMMIE","KENISHA","KAZUKO","KASSIE","GARY","EARLEAN","DAPHINE","CORLISS","CLOTILDE","CAROLYNE","BERNETTA","AUGUSTINA","AUDREA","ANNIS","ANNABELL","YAN","TENNILLE","TAMICA","SELENE","SEAN","ROSANA","REGENIA","QIANA","MARKITA","MACY","LEEANNE","LAURINE","KYM","JESSENIA","JANITA","GEORGINE","GENIE","EMIKO","ELVIE","DEANDRA","DAGMAR","CORIE","COLLEN","CHERISH","ROMAINE","PORSHA","PEARLENE","MICHELINE","MERNA","MARGORIE","MARGARETTA","LORE","KENNETH","JENINE","HERMINA","FREDERICKA","ELKE","DRUSILLA","DORATHY","DIONE","DESIRE","CELENA","BRIGIDA","ANGELES","ALLEGRA","THEO","TAMEKIA","SYNTHIA","STEPHEN","SOOK","SLYVIA","ROSANN","REATHA","RAYE","MARQUETTA","MARGART","LING","LAYLA","KYMBERLY","KIANA","KAYLEEN","KATLYN","KARMEN","JOELLA","IRINA","EMELDA","ELENI","DETRA","CLEMMIE","CHERYLL","CHANTELL","CATHEY","ARNITA","ARLA","ANGLE","ANGELIC","ALYSE","ZOFIA","THOMASINE","TENNIE","SON","SHERLY","SHERLEY","SHARYL","REMEDIOS","PETRINA","NICKOLE","MYUNG","MYRLE","MOZELLA","LOUANNE","LISHA","LATIA","LANE","KRYSTA","JULIENNE","JOEL","JEANENE","JACQUALINE","ISAURA","GWENDA","EARLEEN","DONALD","CLEOPATRA","CARLIE","AUDIE","ANTONIETTA","ALISE","ALEX","VERDELL","VAL","TYLER","TOMOKO","THAO","TALISHA","STEVEN","SO","SHEMIKA","SHAUN","SCARLET","SAVANNA","SANTINA","ROSIA","RAEANN","ODILIA","NANA","MINNA","MAGAN","LYNELLE","LE","KARMA","JOEANN","IVANA","INELL","ILANA","HYE","HONEY","HEE","GUDRUN","FRANK","DREAMA","CRISSY","CHANTE","CARMELINA","ARVILLA","ARTHUR","ANNAMAE","ALVERA","ALEIDA","AARON","YEE","YANIRA","VANDA","TIANNA","TAM","STEFANIA","SHIRA","PERRY","NICOL","NANCIE","MONSERRATE","MINH","MELYNDA","MELANY","MATTHEW","LOVELLA","LAURE","KIRBY","KACY","JACQUELYNN","HYON","GERTHA","FRANCISCO","ELIANA","CHRISTENA","CHRISTEEN","CHARISE","CATERINA","CARLEY","CANDYCE","ARLENA","AMMIE","YANG","WILLETTE","VANITA","TUYET","TINY","SYREETA","SILVA","SCOTT","RONALD","PENNEY","NYLA","MICHAL","MAURICE","MARYAM","MARYA","MAGEN","LUDIE","LOMA","LIVIA","LANELL","KIMBERLIE","JULEE","DONETTA","DIEDRA","DENISHA","DEANE","DAWNE","CLARINE","CHERRYL","BRONWYN","BRANDON","ALLA","VALERY","TONDA","SUEANN","SORAYA","SHOSHANA","SHELA","SHARLEEN","SHANELLE","NERISSA","MICHEAL","MERIDITH","MELLIE","MAYE","MAPLE","MAGARET","LUIS","LILI","LEONILA","LEONIE","LEEANNA","LAVONIA","LAVERA","KRISTEL","KATHEY","KATHE","JUSTIN","JULIAN","JIMMY","JANN","ILDA","HILDRED","HILDEGARDE","GENIA","FUMIKO","EVELIN","ERMELINDA","ELLY","DUNG","DOLORIS","DIONNA","DANAE","BERNEICE","ANNICE","ALIX","VERENA","VERDIE","TRISTAN","SHAWNNA","SHAWANA","SHAUNNA","ROZELLA","RANDEE","RANAE","MILAGRO","LYNELL","LUISE","LOUIE","LOIDA","LISBETH","KARLEEN","JUNITA","JONA","ISIS","HYACINTH","HEDY","GWENN","ETHELENE","ERLINE","EDWARD","DONYA","DOMONIQUE","DELICIA","DANNETTE","CICELY","BRANDA","BLYTHE","BETHANN","ASHLYN","ANNALEE","ALLINE","YUKO","VELLA","TRANG","TOWANDA","TESHA","SHERLYN","NARCISA","MIGUELINA","MERI","MAYBELL","MARLANA","MARGUERITA","MADLYN","LUNA","LORY","LORIANN","LIBERTY","LEONORE","LEIGHANN","LAURICE","LATESHA","LARONDA","KATRICE","KASIE","KARL","KALEY","JADWIGA","GLENNIE","GEARLDINE","FRANCINA","EPIFANIA","DYAN","DORIE","DIEDRE","DENESE","DEMETRICE","DELENA","DARBY","CRISTIE","CLEORA","CATARINA","CARISA","BERNIE","BARBERA","ALMETA","TRULA","TEREASA","SOLANGE","SHEILAH","SHAVONNE","SANORA","ROCHELL","MATHILDE","MARGARETA","MAIA","LYNSEY","LAWANNA","LAUNA","KENA","KEENA","KATIA","JAMEY","GLYNDA","GAYLENE","ELVINA","ELANOR","DANUTA","DANIKA","CRISTEN","CORDIE","COLETTA","CLARITA","CARMON","BRYNN","AZUCENA","AUNDREA","ANGELE","YI","WALTER","VERLIE","VERLENE","TAMESHA","SILVANA","SEBRINA","SAMIRA","REDA","RAYLENE","PENNI","PANDORA","NORAH","NOMA","MIREILLE","MELISSIA","MARYALICE","LARAINE","KIMBERY","KARYL","KARINE","KAM","JOLANDA","JOHANA","JESUSA","JALEESA","JAE","JACQUELYNE","IRISH","ILUMINADA","HILARIA","HANH","GENNIE","FRANCIE","FLORETTA","EXIE","EDDA","DREMA","DELPHA","BEV","BARBAR","ASSUNTA","ARDELL","ANNALISA","ALISIA","YUKIKO","YOLANDO","WONDA","WEI","WALTRAUD","VETA","TEQUILA","TEMEKA","TAMEIKA","SHIRLEEN","SHENITA","PIEDAD","OZELLA","MIRTHA","MARILU","KIMIKO","JULIANE","JENICE","JEN","JANAY","JACQUILINE","HILDE","FE","FAE","EVAN","EUGENE","ELOIS","ECHO","DEVORAH","CHAU","BRINDA","BETSEY","ARMINDA","ARACELIS","APRYL","ANNETT","ALISHIA","VEOLA","USHA","TOSHIKO","THEOLA","TASHIA","TALITHA","SHERY","RUDY","RENETTA","REIKO","RASHEEDA","OMEGA","OBDULIA","MIKA","MELAINE","MEGGAN","MARTIN","MARLEN","MARGET","MARCELINE","MANA","MAGDALEN","LIBRADA","LEZLIE","LEXIE","LATASHIA","LASANDRA","KELLE","ISIDRA","ISA","INOCENCIA","GWYN","FRANCOISE","ERMINIA","ERINN","DIMPLE","DEVORA","CRISELDA","ARMANDA","ARIE","ARIANE","ANGELO","ANGELENA","ALLEN","ALIZA","ADRIENE","ADALINE","XOCHITL","TWANNA","TRAN","TOMIKO","TAMISHA","TAISHA","SUSY","SIU","RUTHA","ROXY","RHONA","RAYMOND","OTHA","NORIKO","NATASHIA","MERRIE","MELVIN","MARINDA","MARIKO","MARGERT","LORIS","LIZZETTE","LEISHA","KAILA","KA","JOANNIE","JERRICA","JENE","JANNET","JANEE","JACINDA","HERTA","ELENORE","DORETTA","DELAINE","DANIELL","CLAUDIE","CHINA","BRITTA","APOLONIA","AMBERLY","ALEASE","YURI","YUK","WEN","WANETA","UTE","TOMI","SHARRI","SANDIE","ROSELLE","REYNALDA","RAGUEL","PHYLICIA","PATRIA","OLIMPIA","ODELIA","MITZIE","MITCHELL","MISS","MINDA","MIGNON","MICA","MENDY","MARIVEL","MAILE","LYNETTA","LAVETTE","LAURYN","LATRISHA","LAKIESHA","KIERSTEN","KARY","JOSPHINE","JOLYN","JETTA","JANISE","JACQUIE","IVELISSE","GLYNIS","GIANNA","GAYNELLE","EMERALD","DEMETRIUS","DANYELL","DANILLE","DACIA","CORALEE","CHER","CEOLA","BRETT","BELL","ARIANNE","ALESHIA","YUNG","WILLIEMAE","TROY","TRINH","THORA","TAI","SVETLANA","SHERIKA","SHEMEKA","SHAUNDA","ROSELINE","RICKI","MELDA","MALLIE","LAVONNA","LATINA","LARRY","LAQUANDA","LALA","LACHELLE","KLARA","KANDIS","JOHNA","JEANMARIE","JAYE","HANG","GRAYCE","GERTUDE","EMERITA","EBONIE","CLORINDA","CHING","CHERY","CAROLA","BREANN","BLOSSOM","BERNARDINE","BECKI","ARLETHA","ARGELIA","ARA","ALITA","YULANDA","YON","YESSENIA","TOBI","TASIA","SYLVIE","SHIRL","SHIRELY","SHERIDAN","SHELLA","SHANTELLE","SACHA","ROYCE","REBECKA","REAGAN","PROVIDENCIA","PAULENE","MISHA","MIKI","MARLINE","MARICA","LORITA","LATOYIA","LASONYA","KERSTIN","KENDA","KEITHA","KATHRIN","JAYMIE","JACK","GRICELDA","GINETTE","ERYN","ELINA","ELFRIEDA","DANYEL","CHEREE","CHANELLE","BARRIE","AVERY","AURORE","ANNAMARIA","ALLEEN","AILENE","AIDE","YASMINE","VASHTI","VALENTINE","TREASA","TORY","TIFFANEY","SHERYLL","SHARIE","SHANAE","SAU","RAISA","PA","NEDA","MITSUKO","MIRELLA","MILDA","MARYANNA","MARAGRET","MABELLE","LUETTA","LORINA","LETISHA","LATARSHA","LANELLE","LAJUANA","KRISSY","KARLY","KARENA","JON","JESSIKA","JERICA","JEANELLE","JANUARY","JALISA","JACELYN","IZOLA","IVEY","GREGORY","EUNA","ETHA","DREW","DOMITILA","DOMINICA","DAINA","CREOLA","CARLI","CAMIE","BUNNY","BRITTNY","ASHANTI","ANISHA","ALEEN","ADAH","YASUKO","WINTER","VIKI","VALRIE","TONA","TINISHA","THI","TERISA","TATUM","TANEKA","SIMONNE","SHALANDA","SERITA","RESSIE","REFUGIA","PAZ","OLENE","NA","MERRILL","MARGHERITA","MANDIE","MAN","MAIRE","LYNDIA","LUCI","LORRIANE","LORETA","LEONIA","LAVONA","LASHAWNDA","LAKIA","KYOKO","KRYSTINA","KRYSTEN","KENIA","KELSI","JUDE","JEANICE","ISOBEL","GEORGIANN","GENNY","FELICIDAD","EILENE","DEON","DELOISE","DEEDEE","DANNIE","CONCEPTION","CLORA","CHERILYN","CHANG","CALANDRA","BERRY","ARMANDINA","ANISA","ULA","TIMOTHY","TIERA","THERESSA","STEPHANIA","SIMA","SHYLA","SHONTA","SHERA","SHAQUITA","SHALA","SAMMY","ROSSANA","NOHEMI","NERY","MORIAH","MELITA","MELIDA","MELANI","MARYLYNN","MARISHA","MARIETTE","MALORIE","MADELENE","LUDIVINA","LORIA","LORETTE","LORALEE","LIANNE","LEON","LAVENIA","LAURINDA","LASHON","KIT","KIMI","KEILA","KATELYNN","KAI","JONE","JOANE","JI","JAYNA","JANELLA","JA","HUE","HERTHA","FRANCENE","ELINORE","DESPINA","DELSIE","DEEDRA","CLEMENCIA","CARRY","CAROLIN","CARLOS","BULAH","BRITTANIE","BOK","BLONDELL","BIBI","BEAULAH","BEATA","ANNITA","AGRIPINA","VIRGEN","VALENE","UN","TWANDA","TOMMYE","TOI","TARRA","TARI","TAMMERA","SHAKIA","SADYE","RUTHANNE","ROCHEL","RIVKA","PURA","NENITA","NATISHA","MING","MERRILEE","MELODEE","MARVIS","LUCILLA","LEENA","LAVETA","LARITA","LANIE","KEREN","ILEEN","GEORGEANN","GENNA","GENESIS","FRIDA","EWA","EUFEMIA","EMELY","ELA","EDYTH","DEONNA","DEADRA","DARLENA","CHANELL","CHAN","CATHERN","CASSONDRA","CASSAUNDRA","BERNARDA","BERNA","ARLINDA","ANAMARIA","ALBERT","WESLEY","VERTIE","VALERI","TORRI","TATYANA","STASIA","SHERISE","SHERILL","SEASON","SCOTTIE","SANDA","RUTHE","ROSY","ROBERTO","ROBBI","RANEE","QUYEN","PEARLY","PALMIRA","ONITA","NISHA","NIESHA","NIDA","NEVADA","NAM","MERLYN","MAYOLA","MARYLOUISE","MARYLAND","MARX","MARTH","MARGENE","MADELAINE","LONDA","LEONTINE","LEOMA","LEIA","LAWRENCE","LAURALEE","LANORA","LAKITA","KIYOKO","KETURAH","KATELIN","KAREEN","JONIE","JOHNETTE","JENEE","JEANETT","IZETTA","HIEDI","HEIKE","HASSIE","HAROLD","GIUSEPPINA","GEORGANN","FIDELA","FERNANDE","ELWANDA","ELLAMAE","ELIZ","DUSTI","DOTTY","CYNDY","CORALIE","CELESTA","ARGENTINA","ALVERTA","XENIA","WAVA","VANETTA","TORRIE","TASHINA","TANDY","TAMBRA","TAMA","STEPANIE","SHILA","SHAUNTA","SHARAN","SHANIQUA","SHAE","SETSUKO","SERAFINA","SANDEE","ROSAMARIA","PRISCILA","OLINDA","NADENE","MUOI","MICHELINA","MERCEDEZ","MARYROSE","MARIN","MARCENE","MAO","MAGALI","MAFALDA","LOGAN","LINN","LANNIE","KAYCE","KAROLINE","KAMILAH","KAMALA","JUSTA","JOLINE","JENNINE","JACQUETTA","IRAIDA","GERALD","GEORGEANNA","FRANCHESCA","FAIRY","EMELINE","ELANE","EHTEL","EARLIE","DULCIE","DALENE","CRIS","CLASSIE","CHERE","CHARIS","CAROYLN","CARMINA","CARITA","BRIAN","BETHANIE","AYAKO","ARICA","AN","ALYSA","ALESSANDRA","AKILAH","ADRIEN","ZETTA","YOULANDA","YELENA","YAHAIRA","XUAN","WENDOLYN","VICTOR","TIJUANA","TERRELL","TERINA","TERESIA","SUZI","SUNDAY","SHERELL","SHAVONDA","SHAUNTE","SHARDA","SHAKITA","SENA","RYANN","RUBI","RIVA","REGINIA","REA","RACHAL","PARTHENIA","PAMULA","MONNIE","MONET","MICHAELE","MELIA","MARINE","MALKA","MAISHA","LISANDRA","LEO","LEKISHA","LEAN","LAURENCE","LAKENDRA","KRYSTIN","KORTNEY","KIZZIE","KITTIE","KERA","KENDAL","KEMBERLY","KANISHA","JULENE","JULE","JOSHUA","JOHANNE","JEFFREY","JAMEE","HAN","HALLEY","GIDGET","GALINA","FREDRICKA","FLETA","FATIMAH","EUSEBIA","ELZA","ELEONORE","DORTHEY","DORIA","DONELLA","DINORAH","DELORSE","CLARETHA","CHRISTINIA","CHARLYN","BONG","BELKIS","AZZIE","ANDERA","AIKO","ADENA","YER","YAJAIRA","WAN","VANIA","ULRIKE","TOSHIA","TIFANY","STEFANY","SHIZUE","SHENIKA","SHAWANNA","SHAROLYN","SHARILYN","SHAQUANA","SHANTAY","SEE","ROZANNE","ROSELEE","RICKIE","REMONA","REANNA","RAELENE","QUINN","PHUNG","PETRONILA","NATACHA","NANCEY","MYRL","MIYOKO","MIESHA","MERIDETH","MARVELLA","MARQUITTA","MARHTA","MARCHELLE","LIZETH","LIBBIE","LAHOMA","LADAWN","KINA","KATHELEEN","KATHARYN","KARISA","KALEIGH","JUNIE","JULIEANN","JOHNSIE","JANEAN","JAIMEE","JACKQUELINE","HISAKO","HERMA","HELAINE","GWYNETH","GLENN","GITA","EUSTOLIA","EMELINA","ELIN","EDRIS","DONNETTE","DONNETTA","DIERDRE","DENAE","DARCEL","CLAUDE","CLARISA","CINDERELLA","CHIA","CHARLESETTA","CHARITA","CELSA","CASSY","CASSI","CARLEE","BRUNA","BRITTANEY","BRANDE","BILLI","BAO","ANTONETTA","ANGLA","ANGELYN","ANALISA","ALANE","WENONA","WENDIE","VERONIQUE","VANNESA","TOBIE","TEMPIE","SUMIKO","SULEMA","SPARKLE","SOMER","SHEBA","SHAYNE","SHARICE","SHANEL","SHALON","SAGE","ROY","ROSIO","ROSELIA","RENAY","REMA","REENA","PORSCHE","PING","PEG","OZIE","ORETHA","ORALEE","ODA","NU","NGAN","NAKESHA","MILLY","MARYBELLE","MARLIN","MARIS","MARGRETT","MARAGARET","MANIE","LURLENE","LILLIA","LIESELOTTE","LAVELLE","LASHAUNDA","LAKEESHA","KEITH","KAYCEE","KALYN","JOYA","JOETTE","JENAE","JANIECE","ILLA","GRISEL","GLAYDS","GENEVIE","GALA","FREDDA","FRED","ELMER","ELEONOR","DEBERA","DEANDREA","DAN","CORRINNE","CORDIA","CONTESSA","COLENE","CLEOTILDE","CHARLOTT","CHANTAY","CECILLE","BEATRIS","AZALEE","ARLEAN","ARDATH","ANJELICA","ANJA","ALFREDIA","ALEISHA","ADAM","ZADA","YUONNE","XIAO","WILLODEAN","WHITLEY","VENNIE","VANNA","TYISHA","TOVA","TORIE","TONISHA","TILDA","TIEN","TEMPLE","SIRENA","SHERRIL","SHANTI","SHAN","SENAIDA","SAMELLA","ROBBYN","RENDA","REITA","PHEBE","PAULITA","NOBUKO","NGUYET","NEOMI","MOON","MIKAELA","MELANIA","MAXIMINA","MARG","MAISIE","LYNNA","LILLI","LAYNE","LASHAUN","LAKENYA","LAEL","KIRSTIE","KATHLINE","KASHA","KARLYN","KARIMA","JOVAN","JOSEFINE","JENNELL","JACQUI","JACKELYN","HYO","HIEN","GRAZYNA","FLORRIE","FLORIA","ELEONORA","DWANA","DORLA","DONG","DELMY","DEJA","DEDE","DANN","CRYSTA","CLELIA","CLARIS","CLARENCE","CHIEKO","CHERLYN","CHERELLE","CHARMAIN","CHARA","CAMMY","BEE","ARNETTE","ARDELLE","ANNIKA","AMIEE","AMEE","ALLENA","YVONE","YUKI","YOSHIE","YEVETTE","YAEL","WILLETTA","VONCILE","VENETTA","TULA","TONETTE","TIMIKA","TEMIKA","TELMA","TEISHA","TAREN","TA","STACEE","SHIN","SHAWNTA","SATURNINA","RICARDA","POK","PASTY","ONIE","NUBIA","MORA","MIKE","MARIELLE","MARIELLA","MARIANELA","MARDELL","MANY","LUANNA","LOISE","LISABETH","LINDSY","LILLIANA","LILLIAM","LELAH","LEIGHA","LEANORA","LANG","KRISTEEN","KHALILAH","KEELEY","KANDRA","JUNKO","JOAQUINA","JERLENE","JANI","JAMIKA","JAME","HSIU","HERMILA","GOLDEN","GENEVIVE","EVIA","EUGENA","EMMALINE","ELFREDA","ELENE","DONETTE","DELCIE","DEEANNA","DARCEY","CUC","CLARINDA","CIRA","CHAE","CELINDA","CATHERYN","CATHERIN","CASIMIRA","CARMELIA","CAMELLIA","BREANA","BOBETTE","BERNARDINA","BEBE","BASILIA","ARLYNE","AMAL","ALAYNA","ZONIA","ZENIA","YURIKO","YAEKO","WYNELL","WILLOW","WILLENA","VERNIA","TU","TRAVIS","TORA","TERRILYN","TERICA","TENESHA","TAWNA","TAJUANA","TAINA","STEPHNIE","SONA","SOL","SINA","SHONDRA","SHIZUKO","SHERLENE","SHERICE","SHARIKA","ROSSIE","ROSENA","RORY","RIMA","RIA","RHEBA","RENNA","PETER","NATALYA","NANCEE","MELODI","MEDA","MAXIMA","MATHA","MARKETTA","MARICRUZ","MARCELENE","MALVINA","LUBA","LOUETTA","LEIDA","LECIA","LAURAN","LASHAWNA","LAINE","KHADIJAH","KATERINE","KASI","KALLIE","JULIETTA","JESUSITA","JESTINE","JESSIA","JEREMY","JEFFIE","JANYCE","ISADORA","GEORGIANNE","FIDELIA","EVITA","EURA","EULAH","ESTEFANA","ELSY","ELIZABET","ELADIA","DODIE","DION","DIA","DENISSE","DELORAS","DELILA","DAYSI","DAKOTA","CURTIS","CRYSTLE","CONCHA","COLBY","CLARETTA","CHU","CHRISTIA","CHARLSIE","CHARLENA","CARYLON","BETTYANN","ASLEY","ASHLEA","AMIRA","AI","AGUEDA","AGNUS","YUETTE","VINITA","VICTORINA","TYNISHA","TREENA","TOCCARA","TISH","THOMASENA","TEGAN","SOILA","SHILOH","SHENNA","SHARMAINE","SHANTAE","SHANDI","SEPTEMBER","SARAN","SARAI","SANA","SAMUEL","SALLEY","ROSETTE","ROLANDE","REGINE","OTELIA","OSCAR","OLEVIA","NICHOLLE","NECOLE","NAIDA","MYRTA","MYESHA","MITSUE","MINTA","MERTIE","MARGY","MAHALIA","MADALENE","LOVE","LOURA","LOREAN","LEWIS","LESHA","LEONIDA","LENITA","LAVONE","LASHELL","LASHANDRA","LAMONICA","KIMBRA","KATHERINA","KARRY","KANESHA","JULIO","JONG","JENEVA","JAQUELYN","HWA","GILMA","GHISLAINE","GERTRUDIS","FRANSISCA","FERMINA","ETTIE","ETSUKO","ELLIS","ELLAN","ELIDIA","EDRA","DORETHEA","DOREATHA","DENYSE","DENNY","DEETTA","DAINE","CYRSTAL","CORRIN","CAYLA","CARLITA","CAMILA","BURMA","BULA","BUENA","BLAKE","BARABARA","AVRIL","AUSTIN","ALAINE","ZANA","WILHEMINA","WANETTA","VIRGIL","VI","VERONIKA","VERNON","VERLINE","VASILIKI","TONITA","TISA","TEOFILA","TAYNA","TAUNYA","TANDRA","TAKAKO","SUNNI","SUANNE","SIXTA","SHARELL","SEEMA","RUSSELL","ROSENDA","ROBENA","RAYMONDE","PEI","PAMILA","OZELL","NEIDA","NEELY","MISTIE","MICHA","MERISSA","MAURITA","MARYLN","MARYETTA","MARSHALL","MARCELL","MALENA","MAKEDA","MADDIE","LOVETTA","LOURIE","LORRINE","LORILEE","LESTER","LAURENA","LASHAY","LARRAINE","LAREE","LACRESHA","KRISTLE","KRISHNA","KEVA","KEIRA","KAROLE","JOIE","JINNY","JEANNETTA","JAMA","HEIDY","GILBERTE","GEMA","FAVIOLA","EVELYNN","ENDA","ELLI","ELLENA","DIVINA","DAGNY","COLLENE","CODI","CINDIE","CHASSIDY","CHASIDY","CATRICE","CATHERINA","CASSEY","CAROLL","CARLENA","CANDRA","CALISTA","BRYANNA","BRITTENY","BEULA","BARI","AUDRIE","AUDRIA","ARDELIA","ANNELLE","ANGILA","ALONA","ALLYN","DOUGLAS","ROGER","JONATHAN","RALPH","NICHOLAS","BENJAMIN","BRUCE","HARRY","WAYNE","STEVE","HOWARD","ERNEST","PHILLIP","TODD","CRAIG","ALAN","PHILIP","EARL","DANNY","BRYAN","STANLEY","LEONARD","NATHAN","MANUEL","RODNEY","MARVIN","VINCENT","JEFFERY","JEFF","CHAD","JACOB","ALFRED","BRADLEY","HERBERT","FREDERICK","EDWIN","DON","RICKY","RANDALL","BARRY","BERNARD","LEROY","MARCUS","THEODORE","CLIFFORD","MIGUEL","JIM","TOM","CALVIN","BILL","LLOYD","DEREK","WARREN","DARRELL","JEROME","FLOYD","ALVIN","TIM","GORDON","GREG","JORGE","DUSTIN","PEDRO","DERRICK","ZACHARY","HERMAN","GLEN","HECTOR","RICARDO","RICK","BRENT","RAMON","GILBERT","MARC","REGINALD","RUBEN","NATHANIEL","RAFAEL","EDGAR","MILTON","RAUL","BEN","CHESTER","DUANE","FRANKLIN","BRAD","RON","ROLAND","ARNOLD","HARVEY","JARED","ERIK","DARRYL","NEIL","JAVIER","FERNANDO","CLINTON","TED","MATHEW","TYRONE","DARREN","LANCE","KURT","ALLAN","NELSON","GUY","CLAYTON","HUGH","MAX","DWAYNE","DWIGHT","ARMANDO","FELIX","EVERETT","IAN","WALLACE","KEN","BOB","ALFREDO","ALBERTO","DAVE","IVAN","BYRON","ISAAC","MORRIS","CLIFTON","WILLARD","ROSS","ANDY","SALVADOR","KIRK","SERGIO","SETH","KENT","TERRANCE","EDUARDO","TERRENCE","ENRIQUE","WADE","STUART","FREDRICK","ARTURO","ALEJANDRO","NICK","LUTHER","WENDELL","JEREMIAH","JULIUS","OTIS","TREVOR","OLIVER","LUKE","HOMER","GERARD","DOUG","KENNY","HUBERT","LYLE","MATT","ALFONSO","ORLANDO","REX","CARLTON","ERNESTO","NEAL","PABLO","LORENZO","OMAR","WILBUR","GRANT","HORACE","RODERICK","ABRAHAM","WILLIS","RICKEY","ANDRES","CESAR","JOHNATHAN","MALCOLM","RUDOLPH","DAMON","KELVIN","PRESTON","ALTON","ARCHIE","MARCO","WM","PETE","RANDOLPH","GARRY","GEOFFREY","JONATHON","FELIPE","GERARDO","ED","DOMINIC","DELBERT","COLIN","GUILLERMO","EARNEST","LUCAS","BENNY","SPENCER","RODOLFO","MYRON","EDMUND","GARRETT","SALVATORE","CEDRIC","LOWELL","GREGG","SHERMAN","WILSON","SYLVESTER","ROOSEVELT","ISRAEL","JERMAINE","FORREST","WILBERT","LELAND","SIMON","CLARK","IRVING","BRYANT","OWEN","RUFUS","WOODROW","KRISTOPHER","MACK","LEVI","MARCOS","GUSTAVO","JAKE","LIONEL","GILBERTO","CLINT","NICOLAS","ISMAEL","ORVILLE","ERVIN","DEWEY","AL","WILFRED","JOSH","HUGO","IGNACIO","CALEB","TOMAS","SHELDON","ERICK","STEWART","DOYLE","DARREL","ROGELIO","TERENCE","SANTIAGO","ALONZO","ELIAS","BERT","ELBERT","RAMIRO","CONRAD","NOAH","GRADY","PHIL","CORNELIUS","LAMAR","ROLANDO","CLAY","PERCY","DEXTER","BRADFORD","DARIN","AMOS","MOSES","IRVIN","SAUL","ROMAN","RANDAL","TIMMY","DARRIN","WINSTON","BRENDAN","ABEL","DOMINICK","BOYD","EMILIO","ELIJAH","DOMINGO","EMMETT","MARLON","EMANUEL","JERALD","EDMOND","EMIL","DEWAYNE","WILL","OTTO","TEDDY","REYNALDO","BRET","JESS","TRENT","HUMBERTO","EMMANUEL","STEPHAN","VICENTE","LAMONT","GARLAND","MILES","EFRAIN","HEATH","RODGER","HARLEY","ETHAN","ELDON","ROCKY","PIERRE","JUNIOR","FREDDY","ELI","BRYCE","ANTOINE","STERLING","CHASE","GROVER","ELTON","CLEVELAND","DYLAN","CHUCK","DAMIAN","REUBEN","STAN","AUGUST","LEONARDO","JASPER","RUSSEL","ERWIN","BENITO","HANS","MONTE","BLAINE","ERNIE","CURT","QUENTIN","AGUSTIN","MURRAY","JAMAL","ADOLFO","HARRISON","TYSON","BURTON","BRADY","ELLIOTT","WILFREDO","BART","JARROD","VANCE","DENIS","DAMIEN","JOAQUIN","HARLAN","DESMOND","ELLIOT","DARWIN","GREGORIO","BUDDY","XAVIER","KERMIT","ROSCOE","ESTEBAN","ANTON","SOLOMON","SCOTTY","NORBERT","ELVIN","WILLIAMS","NOLAN","ROD","QUINTON","HAL","BRAIN","ROB","ELWOOD","KENDRICK","DARIUS","MOISES","FIDEL","THADDEUS","CLIFF","MARCEL","JACKSON","RAPHAEL","BRYON","ARMAND","ALVARO","JEFFRY","DANE","JOESPH","THURMAN","NED","RUSTY","MONTY","FABIAN","REGGIE","MASON","GRAHAM","ISAIAH","VAUGHN","GUS","LOYD","DIEGO","ADOLPH","NORRIS","MILLARD","ROCCO","GONZALO","DERICK","RODRIGO","WILEY","RIGOBERTO","ALPHONSO","TY","NOE","VERN","REED","JEFFERSON","ELVIS","BERNARDO","MAURICIO","HIRAM","DONOVAN","BASIL","RILEY","NICKOLAS","MAYNARD","SCOT","VINCE","QUINCY","EDDY","SEBASTIAN","FEDERICO","ULYSSES","HERIBERTO","DONNELL","COLE","DAVIS","GAVIN","EMERY","WARD","ROMEO","JAYSON","DANTE","CLEMENT","COY","MAXWELL","JARVIS","BRUNO","ISSAC","DUDLEY","BROCK","SANFORD","CARMELO","BARNEY","NESTOR","STEFAN","DONNY","ART","LINWOOD","BEAU","WELDON","GALEN","ISIDRO","TRUMAN","DELMAR","JOHNATHON","SILAS","FREDERIC","DICK","IRWIN","MERLIN","CHARLEY","MARCELINO","HARRIS","CARLO","TRENTON","KURTIS","HUNTER","AURELIO","WINFRED","VITO","COLLIN","DENVER","CARTER","LEONEL","EMORY","PASQUALE","MOHAMMAD","MARIANO","DANIAL","LANDON","DIRK","BRANDEN","ADAN","BUFORD","GERMAN","WILMER","EMERSON","ZACHERY","FLETCHER","JACQUES","ERROL","DALTON","MONROE","JOSUE","EDWARDO","BOOKER","WILFORD","SONNY","SHELTON","CARSON","THERON","RAYMUNDO","DAREN","HOUSTON","ROBBY","LINCOLN","GENARO","BENNETT","OCTAVIO","CORNELL","HUNG","ARRON","ANTONY","HERSCHEL","GIOVANNI","GARTH","CYRUS","CYRIL","RONNY","LON","FREEMAN","DUNCAN","KENNITH","CARMINE","ERICH","CHADWICK","WILBURN","RUSS","REID","MYLES","ANDERSON","MORTON","JONAS","FOREST","MITCHEL","MERVIN","ZANE","RICH","JAMEL","LAZARO","ALPHONSE","RANDELL","MAJOR","JARRETT","BROOKS","ABDUL","LUCIANO","SEYMOUR","EUGENIO","MOHAMMED","VALENTIN","CHANCE","ARNULFO","LUCIEN","FERDINAND","THAD","EZRA","ALDO","RUBIN","ROYAL","MITCH","EARLE","ABE","WYATT","MARQUIS","LANNY","KAREEM","JAMAR","BORIS","ISIAH","EMILE","ELMO","ARON","LEOPOLDO","EVERETTE","JOSEF","ELOY","RODRICK","REINALDO","LUCIO","JERROD","WESTON","HERSHEL","BARTON","PARKER","LEMUEL","BURT","JULES","GIL","ELISEO","AHMAD","NIGEL","EFREN","ANTWAN","ALDEN","MARGARITO","COLEMAN","DINO","OSVALDO","LES","DEANDRE","NORMAND","KIETH","TREY","NORBERTO","NAPOLEON","JEROLD","FRITZ","ROSENDO","MILFORD","CHRISTOPER","ALFONZO","LYMAN","JOSIAH","BRANT","WILTON","RICO","JAMAAL","DEWITT","BRENTON","OLIN","FOSTER","FAUSTINO","CLAUDIO","JUDSON","GINO","EDGARDO","ALEC","TANNER","JARRED","DONN","TAD","PRINCE","PORFIRIO","ODIS","LENARD","CHAUNCEY","TOD","MEL","MARCELO","KORY","AUGUSTUS","KEVEN","HILARIO","BUD","SAL","ORVAL","MAURO","ZACHARIAH","OLEN","ANIBAL","MILO","JED","DILLON","AMADO","NEWTON","LENNY","RICHIE","HORACIO","BRICE","MOHAMED","DELMER","DARIO","REYES","MAC","JONAH","JERROLD","ROBT","HANK","RUPERT","ROLLAND","KENTON","DAMION","ANTONE","WALDO","FREDRIC","BRADLY","KIP","BURL","WALKER","TYREE","JEFFEREY","AHMED","WILLY","STANFORD","OREN","NOBLE","MOSHE","MIKEL","ENOCH","BRENDON","QUINTIN","JAMISON","FLORENCIO","DARRICK","TOBIAS","HASSAN","GIUSEPPE","DEMARCUS","CLETUS","TYRELL","LYNDON","KEENAN","WERNER","GERALDO","COLUMBUS","CHET","BERTRAM","MARKUS","HUEY","HILTON","DWAIN","DONTE","TYRON","OMER","ISAIAS","HIPOLITO","FERMIN","ADALBERTO","BO","BARRETT","TEODORO","MCKINLEY","MAXIMO","GARFIELD","RALEIGH","LAWERENCE","ABRAM","RASHAD","KING","EMMITT","DARON","SAMUAL","MIQUEL","EUSEBIO","DOMENIC","DARRON","BUSTER","WILBER","RENATO","JC","HOYT","HAYWOOD","EZEKIEL","CHAS","FLORENTINO","ELROY","CLEMENTE","ARDEN","NEVILLE","EDISON","DESHAWN","NATHANIAL","JORDON","DANILO","CLAUD","SHERWOOD","RAYMON","RAYFORD","CRISTOBAL","AMBROSE","TITUS","HYMAN","FELTON","EZEQUIEL","ERASMO","STANTON","LONNY","LEN","IKE","MILAN","LINO","JAROD","HERB","ANDREAS","WALTON","RHETT","PALMER","DOUGLASS","CORDELL","OSWALDO","ELLSWORTH","VIRGILIO","TONEY","NATHANAEL","DEL","BENEDICT","MOSE","JOHNSON","ISREAL","GARRET","FAUSTO","ASA","ARLEN","ZACK","WARNER","MODESTO","FRANCESCO","MANUAL","GAYLORD","GASTON","FILIBERTO","DEANGELO","MICHALE","GRANVILLE","WES","MALIK","ZACKARY","TUAN","ELDRIDGE","CRISTOPHER","CORTEZ","ANTIONE","MALCOM","LONG","KOREY","JOSPEH","COLTON","WAYLON","VON","HOSEA","SHAD","SANTO","RUDOLF","ROLF","REY","RENALDO","MARCELLUS","LUCIUS","KRISTOFER","BOYCE","BENTON","HAYDEN","HARLAND","ARNOLDO","RUEBEN","LEANDRO","KRAIG","JERRELL","JEROMY","HOBERT","CEDRICK","ARLIE","WINFORD","WALLY","LUIGI","KENETH","JACINTO","GRAIG","FRANKLYN","EDMUNDO","SID","PORTER","LEIF","JERAMY","BUCK","WILLIAN","VINCENZO","SHON","LYNWOOD","JERE","HAI","ELDEN","DORSEY","DARELL","BRODERICK","ALONSO"};
    numNames = countArray(names);
    selSort(names, 0, numNames - 1);
    unsigned long sum;
    int i, score;
    sum = 0;
    for(i = 0; i < numNames; i++)
        sum += (nameScore(names[i]) * (i + 1));
    printf("Total Score: %u", sum);
    return 0;
}

int nameScore(char *s)
{
    int i, sum;
    sum = 0;
    for(i = 0; i < strlen(s); i++)
        sum += (s[i] - 'A' + 1);
    return sum;
}

int countArray(char *names[])
{
    int i;
    for(i = 0; strcmp(names[i], "ALONSO"); i++)
    ;
    return i + 1;
}

int selSort(char *names[], int start, int end)
{
    char *sel;
    sel = names[start];
    if (end <= start)
        return 1;
    char *sorted[numNames];
    int i;
    int front, back;
    front = start;
    back = end;
    for(i = 0; i < start; i++)
        sorted[i] = names[i];
    for(i = start + 1; i <= end; i++){
        if (strcmp(names[i], sel) > 0)
            sorted[back--] = names[i];
            else sorted[front++] = names[i];
    }
    sorted[front] = sel;
    for(i = end + 1; i < numNames; i++)
        sorted[i] = names[i];
    for(i = 0; i < numNames; i++)
        names[i] = sorted[i];
    selSort(names, start, front - 1);
    selSort(names, front + 1, end);
    return 1;
}

int printNames(char *names[])
{
    int i;
    for(i = 0; i < numNames; i++){
        printf("%s\t", names[i]);
        if ((i + 1) % 5 == 0)
            printf("\n");
    }
    return 1;
}

/*

PROBLEM #17

#define NUM_WORDS 1000

int main()
{
    int i, sum, ones, tens, huns;
    sum = 0;
    for(i = 1; i < NUM_WORDS; i++){
        ones = i % 10;
        tens = (i / 10) % 10;
        huns = i / 100;
        switch (huns)
        {
            case 0:
                break;
            case 1:
            case 2:
            case 6:
                sum += 10;
                break;
            case 3:
            case 7:
            case 8:
                sum += 12;
                break;
            case 4:
            case 5:
            case 9:
                sum += 11;
                break;
        }
        if(i % 100 == 0)
            continue;
        if(huns > 0)
            sum += 3;
        switch (tens)
        {
            case 0:
                break;
            case 1:
                switch (ones)
                {
                    case 0:
                        sum += 3;
                        break;
                    case 1:
                    case 2:
                        sum += 6;
                        break;
                    case 3:
                    case 4:
                    case 8:
                    case 9:
                        sum += 8;
                        break;
                    case 5:
                    case 6:
                        sum += 7;
                        break;
                    case 7:
                        sum += 9;
                        break;
                }
                break;
            case 2:
            case 3:
            case 8:
            case 9:
                sum += 6;
                break;
            case 4:
            case 5:
            case 6:
                sum += 5;
                break;
            case 7:
                sum += 7;
                break;
        }
        if (tens == 1)
            continue;
        switch (ones)
        {
            case 0:
                break;
            case 1:
            case 2:
            case 6:
                sum += 3;
                break;
            case 3:
            case 7:
            case 8:
                sum += 5;
                break;
            case 4:
            case 5:
            case 9:
                sum += 4;
                break;
        }
    }
    sum += 11;
    printf("Number of letters: %d\n", sum);
    return 0;
}

/*

PROBLEM #21

#define MAX 10000

int main()
{
    int numFact[MAX];
    int found[MAX];
    numFact[0] = 0;
    numFact[1] = 0;
    int i, j;
    for(i = 0; i < MAX; i++)
        found[i] = 0;
    for(i = 2; i < MAX; i++)
        numFact[i] = sumFact(i);
    int sum;
    sum = 0;
    j = 0;
    for(i = 2; i < MAX; i++)
        if(numFact[i] < MAX)
            if(numFact[numFact[i]] == i && i != numFact[i]){
                sum += (i + numFact[i]);
                found[j++] = i;
                found[j++] = numFact[i];
                numFact[i] = 0;
            }
    printf("Amicable Sum: %d\n", sum);
    printFound(found);
}

int printFound(int x[])
{
    int i;
    for(i = 0; i < MAX; i++){
        if (x[i] == 0)
            break;
        printf("%d\n", x[i]);
    }
    return 1;
}

int sumFact(int x)
{
    int i, sum;
    sum = 1;
    for(i = 2; i <= sqrt(x); i++){
        if (x % i == 0){
            sum += i;
            if (i != sqrt(x))
                sum += (x / i);
        }
    }
    return sum;
}

/*

PROBLEM #18 AND #67

struct Node {
   int value;
   int greatestPath;
   int pos;
};

typedef struct Node Node;

Node *initNode(Node *node, int pos);
Node *updatePath(Node *node, Node *nodes[]);

int main()
{
    int rows[5050];
    inputArray(rows);
    Node *nodes[5050];
    Node *cur;
    int i;
    for(i = 0; i < 5050; i++)
        nodes[i] = initNode(nodes[i], i);
    for(i = 0; i < 5050; i++)
        nodes[i]->value = rows[i];
    for(i = 5049; i >= 0; i--)
        nodes[i] = updatePath(nodes[i], nodes);
    printf("Longest Path: %d\n", nodes[0]->greatestPath);
    return 0;
}

int inputArray(int a[])
{
    FILE *fp;
    fp = fopen("triangle.txt", "r");
    int cur, i;
    for(i = 0; i < 5050; i++)
    {
        fscanf(fp, "%d", &cur);
        a[i] = cur;
    }
    fclose(fp);
    return 1;
}

Node *updatePath(Node *node, Node *nodes[])
{
    if (node->pos >= 4950)
    {
        node->greatestPath = node->value;
        return node;
    }
    int row, col, pos, g1, g2;
    pos = node->pos;
    row = findRow(pos);
    col = findCol(pos);
    g1 = nodes[findPos(row + 1, col)]->greatestPath;
    g2 = nodes[findPos(row + 1, col + 1)]->greatestPath;
    if (g1 > g2)
        node->greatestPath = g1 + node->value;
        else node->greatestPath = g2 + node->value;
    return node;
}

Node *initNode(Node *node, int pos)
{
    node = malloc(sizeof(Node));
    node->value = 0;
    node->greatestPath = 0;
    node->pos = pos;
    return node;
}

int findRow(int pos)
{
    if (pos == 0)
        return 0;
    int i;
    for(i = 1; i <= 100; i++)
        if (pos < triangle(i))
            return --i;
}

int findCol(int pos)
{
    return pos - triangle(findRow(pos));
}

int findPos(int row, int col)
{
    return triangle(row) + col;
}

int triangle(int x)
{
    int i, sum;
    sum = 0;
    for(i = 1; i <= x; i++)
        sum += i;
    return sum;
}
/*

int main()
{
    int rows[15][15] =
    {
    {75, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,},
    {95, 64, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,},
    {17, 47, 82, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,},
    {18, 35, 87, 10, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,},
    {20, 4, 82, 47, 65, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,},
    {19, 1, 23, 75, 3, 34, 0,  0,  0,  0,  0,  0,  0,  0,  0,},
    {88, 2, 77, 73, 7, 63, 67, 0,  0,  0,  0,  0,  0,  0,  0,},
    {99, 65, 4, 28, 6, 16, 70, 92, 0,  0,  0,  0,  0,  0,  0,},
    {41, 41, 26, 56, 83, 40, 80, 70, 3,  0,  0,  0,  0,  0,  0,},
    {41, 48, 72, 33, 47, 32, 37, 16, 94, 29, 0,  0,  0,  0,  0,},
    {53, 71, 44, 65, 25, 43, 91, 52, 97, 51, 14, 0,  0,  0,  0,},
    {70, 11, 33, 28, 77, 73, 17, 78, 39, 68, 17, 57, 0,  0,  0,},
    {91, 71, 52, 38, 17, 14, 91, 43, 58, 50, 27, 29, 48, 0,  0,},
    {63, 66, 4, 68, 89, 53, 67, 30, 73, 16, 69, 87, 40, 31, 0,},
    {04, 62, 98, 27, 23, 9, 70, 98, 73, 93, 38, 53, 60, 04, 23}
    };

}

/*

PROBLEM #48

#define END_TERM 1000
#define MAX_DIGITS 10

int main()
{
    int cur[MAX_DIGITS];
    int sum[MAX_DIGITS];
    int i, j;
    for(i = 0; i < MAX_DIGITS; i++)
        sum[i] = 0;
    for(i = 1; i <= END_TERM; i++){
        cur[0] = i;
        for(j = 1; j < MAX_DIGITS; j++)
            cur[j] = 0;
        powArray(cur, i, i);
        add10(cur, sum);
    }
    printArrayNum10(sum);
    return 0;
}

int add10(int a[], int b[])
{
    int i, sum;
    for(i = 0; i < MAX_DIGITS; i++)
    {
        sum = a[i] + b[i];
        b[i] = sum % 10;
        if (i == MAX_DIGITS - 1)
            return 1;
        if (sum > 9)
            b[i + 1] += sum / 10;
    }
    return 1;
}

int powArray(int a[], int pow, int x)
{
    pow--;
    for(pow; pow > 0; pow--)
        mult10(a, x, 0, 0);
    return 1;
}

int mult10(int a[], int x, int pos, int carry)
{
    int prod;
    prod = a[pos] * x;
    prod += carry;
    a[pos] = prod % 10;
    if (pos == MAX_DIGITS - 1)
        return 1;
    mult10(a, x, pos + 1, prod / 10);
}

int printArrayNum10(int x[])
{
    int i;
    for(i = MAX_DIGITS - 1; i >= 0; i--)
        printf("%d", x[i]);
    printf("\n");
    return 1;
}

/*

PROBLEM #25

#define MAX_DIGITS 1002
#define DES_DIGITS 1000

int main()
{
    int a[MAX_DIGITS];
    int b[MAX_DIGITS];
    int i;
    for(i = 0; i < MAX_DIGITS; i++)
        a[i] = b[i] = 0;
    a[0] = b[0] = 1;
    int par;
    par = 0;
    unsigned long count;
    count = 2;
    while(a[DES_DIGITS - 1] == 0 &&
          b[DES_DIGITS - 1] == 0)
    {
        if (par)
            add(a, b);
            else add(b, a);
        count++;
        par = 1 - par;
    }
    printf("Term: %u\n", count);
    printArrayNum(a);
    printArrayNum(b);
    return 0;
}

int add(int a[], int b[])
{
    int i, sum;
    for(i = 0; i < MAX_DIGITS; i++)
    {
        sum = a[i] + b[i];
        if (sum > 9)
            b[i + 1] += sum / 10;
        b[i] = sum % 10;
    }
    return 1;
}

int printArrayNum(int x[])
{
    int i;
    for(i = MAX_DIGITS - 1; x[i] == 0; i--)
    ;
    for(i; i >= 0; i--)
        printf("%d", x[i]);
    printf("\n");
    return 1;
}

/*

PROBLEM #15

long long paths(int row, int col);

#define GRID_SIZE 20

long long pathLengths[GRID_SIZE + 1][GRID_SIZE + 1];

int main()
{
    int i, j;
    for(i = 0; i <= GRID_SIZE; i++)
        for(j = 0; j <= GRID_SIZE; j++)
            pathLengths[i][j] = 0;
    for(i = 0; i <= GRID_SIZE; i++)
        for(j = 0; j <= GRID_SIZE; j++)
            paths(i, j);
    printf("Paths: %lld\n", pathLengths[GRID_SIZE][GRID_SIZE]);
    return 0;
}

long long paths(int row, int col)
{
    if(pathLengths[row][col])
        return pathLengths[row][col];
    if(row == 0 || col == 0)
    {
        pathLengths[row][col] = 1;
        return 1;
    }
    long long length;
    length = paths(row - 1, col) + paths(row, col - 1);
    pathLengths[row][col] = length;
    pathLengths[col][row] = length;
    return length;
}

/*

unsigned long numSplits(int target, int divisions);
int fillStrats(unsigned long strats[]);
unsigned long bruteForce();
unsigned long capMethod(int cap, int depth);

#define GRID_SIZE 20

int main()
{
    int strats[GRID_SIZE + 1];
    strats[0] = 0;
    fillStrats(strats);
    unsigned long sum;
    sum = 0;
    int i;
    for(i = 1; i <= GRID_SIZE; i++)
    {
        sum += strats[i] * strats[i];
        sum += strats[i] * strats[i - 1];
    }
    sum *= 2;
    printf("Number of paths: %u\n", sum);
    printf("Cap method: %u\n", capMethod(GRID_SIZE, 1));
    printf("Brute force: %u\n", bruteForce());
    return 0;
}

unsigned long capMethod(int cap, int depth)
{
    if (depth == GRID_SIZE)
        return cap + 1;
    unsigned long sum;
    sum = 0;
    int i;
    for(i = 0; i <= cap; i++)
        sum += capMethod(i, depth + 1);
    return sum;
}

int fillStrats(unsigned long strats[])
{
    int i;
    for(i = 1; i <= GRID_SIZE; i++)
        strats[i] = numSplits(GRID_SIZE, i);
    return 1;
}

unsigned long numSplits(int target, int divisions)
{
    if(divisions > target || divisions < 1)
        return 0;
    if(divisions == 1)
        return 1;
    int i;
    unsigned long sum;
    sum = 0;
    for(i = 1; i <= target - divisions + 1; i++)
        sum += numSplits(target - i, divisions - 1);
    return sum;
}

int search(int nodes[], int pos);

unsigned long pathCount;

unsigned long bruteForce()
{
    int numNodes;
    numNodes = (GRID_SIZE + 1) * (GRID_SIZE + 1);
    int nodes[numNodes];
    int i;
    for(i = 0; i < numNodes; i++)
        nodes[i] = 0;
    pathCount = 0;
    search(nodes, 0);
    return pathCount;
}

int search(int nodes[], int pos)
{
    int numNodes;
    numNodes = (GRID_SIZE + 1) * (GRID_SIZE + 1);
    if (pos == numNodes - 1)
    {
        pathCount++;
        return 1;
    }
    nodes[pos] = 1;
    int hasR, hasD;
    int r, d;
    r = pos + 1;
    d = pos + (GRID_SIZE + 1);
    hasR = hasD = 0;
    if ((pos + 1) % (GRID_SIZE + 1) != 0)
        hasR = 1;
    if (pos < (GRID_SIZE + 1) * ((GRID_SIZE + 1) - 1))
        hasD = 1;
    if (hasR && !nodes[r])
        search(nodes, r);
    if (hasD && !nodes[d])
        search(nodes, d);
    nodes[pos] = 0;
    return 1;
}

/*

PROBLEM #12

int fillTriangles(unsigned long triangles[]);
int numFactors(unsigned long x);

#define BIG_TRIANGLE 100000

int main()
{
    unsigned long triangles[BIG_TRIANGLE];
    fillTriangles(triangles);
    int i, curFact;
    for(i = 1; i < BIG_TRIANGLE; i++)
    {
        curFact = numFactors(triangles[i]);
        if (curFact >= 500)
            break;
    }
    if(i == BIG_TRIANGLE)
        i--;
    printf("Number: %u\n", triangles[i]);
    printf("Factors: %d\n", curFact);
    return 0;
}

int fillTriangles(unsigned long triangles[])
{
    int i;
    triangles[0] = 0;
    for(i = 1; i < BIG_TRIANGLE; i++)
        triangles[i] = triangles[i - 1] + i;
    return 1;
}

int numFactors(unsigned long x)
{
    int i, factors;
    factors = 0;
    double root;
    root = sqrt(x);
    for(i = 1; i <= root; i++)
        if (x % i == 0)
            factors += 2;
    if (root == round(root))
        factors--;
    return factors;
}



PROBLEM #14

unsigned long chainLength(unsigned long x);
int printChain(unsigned long x);

#define MAX_SIZE 100000000
#define RANGE_MIN 500000
#define RANGE_MAX 1000000

int lengths[MAX_SIZE];

int main(int argc, char *argv[])
{
    unsigned long i, best, cur, bestChain;
    best = bestChain = cur = 0;
    lengths[0] = 0;
    lengths[1] = 1;
    for (i = 2; i < MAX_SIZE; i++)
        lengths[i] = 0;
    for (i = RANGE_MIN; i < RANGE_MAX; i++){
        cur = chainLength(i);
        if (cur > bestChain){
            bestChain = cur;
            best = i;
        }
    }
    printf("%ld\n\n", best);
    printChain(best);
    printf("\nNumber of steps: %d\n", bestChain);
}

unsigned long chainLength(unsigned long x)
{
    if (x < MAX_SIZE)
        if (lengths[x])
            return lengths[x];
    unsigned long nextLength;
    if (x % 2 == 0)
        nextLength = chainLength(x / 2) + 1;
    else nextLength = chainLength(3 * x + 1) + 1;
    if (x < MAX_SIZE)
        lengths[x] = nextLength;
    return nextLength;
}

int printChain(unsigned long x)
{
    int count;
    count = 1;
    while (x != 1)
    {
        printf("%u\t", x);
        if (x % 2 == 0)
            x = x / 2;
        else x = 3 * x + 1;
        if (count % 4 == 0)
            printf("\n");
        count++;
    }
    printf("%u\n\n", x);
    printf("Printed steps: %d\n", count);
    return 1;
}


/*

PROBLEM #20 (Unsuccessful)
int DSIZE;

void mult(int digits[], int x);
void addto(int digits[], int i, int val);

int main(int argc, char *argv[])
{
    int i, fact;
    fact = atoi(argv[1]);
    DSIZE = 2 * fact;
    int sum;
    int digits[DSIZE];
    for (i = 0; i < DSIZE; i++)
        digits[i] = 0;
    digits[0] = 1;
    for (i = 1; i <= fact; i++)
        mult(digits, i);
    sum = 0;
    for (i = 0; i < DSIZE; i++)
        sum += digits[i];
    printf("%d\n", sum);
}

void mult(int digits[], int x)
{
    int i, carry, prod, ones, tens;
    int row1[DSIZE], row2[DSIZE];
    ones = tens = carry = 0;
    if (x < 10)
        ones = x;
    else {
    ones = x % 10;
    tens = x / 10;
    }
    for (i = 0; i < DSIZE; i++)
    {
        prod = (digits[i] * ones) + carry;
        row1[i] = prod % 10;
        carry = prod / 10;
    }
    carry = 0;
    row2[0] = 0;
    for (i = 0; i < DSIZE; i++)
    {
        prod = (digits[i] * tens) + carry;
        row2[i + 1] = prod % 10;
        carry = prod / 10;
    }
    for (i = DSIZE - 1; i >= 0; i--)
        digits[i] = row1[i];
    for (i = DSIZE - 1; i >= 0; i--)
        addto(digits, i, row2[i]);
}

void addto(int digits[], int i, int val)
{
    if (digits[i] + val > 9)
    {
        digits[i] = (digits[i] + val) % 10;
        digits[i + 1]++;
    }
    else digits[i] += val;
}

/*

//PROBLEM #16

#define DSIZE 400
#define POWER 1000

void doub(int digits[]);

int main(int argc, char *argv[])
{
    int i;
    long sum;
    int digits[DSIZE];
    for (i = 0; i < DSIZE; i++)
        digits[i] = 0;
    digits[0] = 2;
    for (i = 1; i < POWER; i++)
        doub(digits);
    sum = 0;
    for (i = 0; i < DSIZE; i++)
        sum += digits[i];
    printf("%ld\n", sum);
}

void doub(int digits[])
{
    int i;
    for (i = DSIZE - 1; i >= 0; i--)
    {
        if (digits[i] * 2 > 9)
        {
            digits[i] = (digits[i] * 2) % 10;
            (digits[i + 1])++;
        }
        else digits[i] *= 2;
    }
}

/*
//PROBLEM 12

int numFactors(unsigned int x);

int main(int argc, char *argv[])
{
    unsigned int i, triangle;
    int factors;
    for (i = 1; ; i++)
    {
        triangle = (i * (i + 1)) / 2;
        factors = numFactors(i);
        if (factors >= atoi(argv[1]))
            break;
    }
    printf("TRIANGLE NUMBER: %d\n", i);
    printf("VALUE: %d\n", triangle);
    printf("NUMBER OF FACTORS: %d\n", factors);
}

int numFactors(unsigned int x)
{
    int prod, i, count;
    unsigned int cur1, cur2, half;
    i = 2;
    if (x % 2 == 0)
    {
        cur1 = x / 2;
        cur2 = x + 1;
    }
    else
    {
        cur1 = x;
        cur2 = (x + 1) / 2;
    }
    half = (cur1 * cur2) / 2;
    prod = 1;
    while (i <= half)
    {
        if (cur1 % i == 0)
        {
            count++;
            cur1 = cur1 / i;
        }
        else if (cur2 % i == 0)
        {
            count++;
            cur2 = cur2 / i;
        }
        else
        {
            prod *= (count + 1);
            count = 0;
            i++;
        }
    }
    if (prod == 1)
        return 2;
    else return prod;
}

*/

