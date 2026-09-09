// Réutilise le noyau d'intervalles figé, sans le modifier.
#include "inertie_binaire_native.cpp"

extern "C" int vinkulum_compare_mass(
    I64 n,const I64* mp,const I64* mj,const double* mv,const double* ell,
    double alpha,double beta,const I64* permutation,I64 budget_ops,I64 budget_coeff,
    Pivot* lower_pivots,Pivot* upper_pivots,Summary* summaries,I64* diagnostic,
    double* bad_interval) noexcept {
    Arithmetic ar(budget_ops);I64 phase=0;
    try {
        if(!environment_valid())throw Failure(8);
        std::vector<I64> order(permutation,permutation+n);
        for(phase=0;phase<2;++phase) {
            Matrix a(n,budget_coeff);
            for(I64 i=0;i<n;++i)for(I64 k=mp[i];k<mp[i+1];++k)if(mj[k]>=i) {
                auto point=Arithmetic::point(mv[k]);
                a.set(i,mj[k],phase==0?point:Arithmetic::neg(point));
            }
            for(I64 i=0;i<n;++i) {
                auto term=ar.mul(Arithmetic::point(phase==0?alpha:beta),Arithmetic::point(ell[i]));
                a.set(i,i,phase==0?ar.sub(a.get(i,i),term):ar.add(a.get(i,i),term));
            }
            eliminate(a,ar,order,true,phase==0?lower_pivots:upper_pivots,summaries[phase]);
        }
        diagnostic[0]=ar.operations;return 0;
    } catch(const Failure& e) {
        diagnostic[0]=ar.operations;diagnostic[1]=phase;diagnostic[2]=e.index;
        bad_interval[0]=e.interval.lo;bad_interval[1]=e.interval.hi;return static_cast<int>(e.code);
    } catch(...) {diagnostic[0]=ar.operations;diagnostic[1]=phase;return 9;}
}
