#ifndef USEFUL_H_INCLUDED
#define USEFUL_H_INCLUDED

int divideAway(int n, int d)
{
    if(n % d == 0)
        return divideAway(n/d, d);
    else return n;
}

int numToArray(int x, int a[], int max)
{
    int i;
    for(i = 0; i < max; i++)
        a[i] = (x / (expt(10, i))) % 10;
    for(i = max - 1; a[i] == 0; i--)
        a[i] = -1;
    return 1;
}

int arrayToNum(int a[], int max)
{
    int i, sum;
    sum = i = 0;
    while(i < max && a[i] != -1)
        sum += (a[i] * expt(10, i++));
    return sum;
}

int printArrayNum(int x[], int max)
{
    int i;
    for(i = max - 1; i >= 0; i--)
        printf("%d", x[i]);
    printf("\n");
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

int expt(int x, int y)
{
    int prod;
    prod = 1;
    for(y; y > 0; y--)
        prod *= x;
    return prod;
}

unsigned long expt_u(int x, int y)
{
    unsigned long prod;
    prod = 1;
    for(y; y > 0; y--)
        prod *= x;
    return prod;
}

int isNear(double a, double b, double width)
{
    double diff;
    diff = a - b;
    if(diff < width && diff > width*-1)
        return 1;
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

#endif // USEFUL_H_INCLUDED
