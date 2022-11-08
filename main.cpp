#include <iostream>
#include <cmath>
#include <cinttypes>

using namespace std;

#define MAXS 30001

int sumpropdiv[MAXS];

int abundant[MAXS];
int nabundant;

uint64_t is_abundant_sum[1 + (MAXS>>6)];

uint64_t iscompound[1 + (MAXS>>6)];

#define baisset(ba, i) ((ba)[(i)>>6] & (1ull << ((i) & 63))) != 0
#define baisclear(ba, i) ((ba)[(i)>>6] & (1ull << ((i) & 63))) == 0
#define baset(ba, i) (ba)[(i)>>6] |= (1ull << ((i) & 63))

void fillprimes(int upperBound) {
	int i, j;

	baset(iscompound, 0);
	baset(iscompound, 1);
	int sqrtUb = sqrt(upperBound);
	for (i = 2; i <= sqrtUb; i++) {
		if (baisclear(iscompound, i)) {
            sumpropdiv[i] = -1;
			for (j = i*i; j < upperBound; j += i) {
				baset(iscompound, j);
                sumpropdiv[j] = i;
			}
		}
	}

	for (i = sqrtUb + 1; i < upperBound; i++) {
		if (baisclear(iscompound, i)) {
            sumpropdiv[i] = -1;
		}
	}
}

void fill_sum_of_prop_div(int limit)
{
    //First: sum of div
    sumpropdiv[1] = 1;
    for(int i = 2; i < limit; i++){
        int p = sumpropdiv[i];
        if(p == -1){
            //Prime
            sumpropdiv[i] = i+1;
        }else{
            int div = p*p;
            int sdiv = 1+p;
            while(i%div == 0){
                sdiv += div;
                div*=p;
            }
            sumpropdiv[i] = sdiv*sumpropdiv[(i*p)/div];
        }
    }

    //sum of prop div: sum of div(n) - n
    for(int i = 1; i < limit; i++){
        sumpropdiv[i] -= i;
    }
}

void fill_abundant(int limit)
{
    for(int i = 1; i < limit; i++)
    {
        if(sumpropdiv[i] > i){
            abundant[nabundant++] = i;
        }
    }
}

int main()
{
    fillprimes(MAXS);
    fill_sum_of_prop_div(MAXS);    
    fill_abundant(MAXS);

    for(int i = 0; i < nabundant; i++){
        for(int j = 0; j < nabundant; j++){
            int absum = abundant[i]+abundant[j];
            if(absum >= MAXS) continue;

            baset(is_abundant_sum, absum);
        }
    }

    int64_t res = 0;

    for(int i = 1; i < MAXS; i++){
        if(baisclear(is_abundant_sum, i)){
            res += i;
        }
    }

    cout << res << endl;

}