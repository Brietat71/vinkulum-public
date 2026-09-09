// Prototype de recherche : congruences intervalles sur les entrées binary64.
// Compiler sans fast-math et sans contraction FMA ; aucun mode global modifié.
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>
#ifdef __FAST_MATH__
#error Fast math invalide le contrat des intervalles
#endif

using I64 = std::int64_t;
static bool environment_valid() {
    static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53);
    volatile double tiny=std::numeric_limits<double>::denorm_min(),normal=std::numeric_limits<double>::min(),half=.5;
    double sum=tiny+tiny,product=normal*half;
    std::uint64_t sum_bits,product_bits;
    std::memcpy(&sum_bits,&sum,sizeof(double));
    std::memcpy(&product_bits,&product,sizeof(double));
    // Comparaison des bits : DAZ pourrait aussi fausser une comparaison FP.
    return std::fegetround()==FE_TONEAREST && sum_bits==2 && product_bits==0x0008000000000000ULL;
}
struct Interval { double lo, hi; };
struct Pivot {
    I64 i, j, positive, negative;
    double alo, ahi, blo, bhi, clo, chi, dlo, dhi;
};
struct Summary {
    I64 positive, negative, zero, count, width, peak, created, initial;
};
struct Failure : std::runtime_error {
    I64 code, index;
    Interval interval;
    Failure(I64 c, I64 i=-1, Interval a={0,0})
        : std::runtime_error("certificat indisponible"), code(c), index(i), interval(a) {}
};
struct Arithmetic {
    I64 operations=0, budget;
    explicit Arithmetic(I64 b): budget(b) {}
    void count(I64 n) {
        if (operations > budget-n) throw Failure(2);
        operations += n;
    }
    static bool zero(Interval a) { return a.lo==0 && a.hi==0; }
    static bool one(Interval a) { return a.lo==1 && a.hi==1; }
    static Interval point(double a) { return {a,a}; }
    static Interval neg(Interval a) { return {-a.hi,-a.lo}; }
    static Interval enclosure(double lo,double hi) {
        lo=std::nextafter(lo,-std::numeric_limits<double>::infinity());
        hi=std::nextafter(hi,std::numeric_limits<double>::infinity());
        if (!std::isfinite(lo)||!std::isfinite(hi)) throw Failure(3);
        return {lo,hi};
    }
    Interval add(Interval a,Interval b) {
        if(zero(a)) return b;
        if(zero(b)) return a;
        count(2);return enclosure(a.lo+b.lo,a.hi+b.hi);
    }
    Interval sub(Interval a,Interval b) {
        if(zero(b)) return a;
        if(zero(a)) return neg(b);
        count(2);return enclosure(a.lo-b.hi,a.hi-b.lo);
    }
    Interval mul(Interval a,Interval b) {
        if(zero(a)||zero(b)) return {0,0};
        if(one(a)) return b;
        if(one(b)) return a;
        double lo,hi;
        if(a.lo>=0) {
            if(b.lo>=0) {lo=a.lo*b.lo;hi=a.hi*b.hi;}
            else if(b.hi<=0) {lo=a.hi*b.lo;hi=a.lo*b.hi;}
            else {lo=a.hi*b.lo;hi=a.hi*b.hi;}
        } else if(a.hi<=0) {
            if(b.lo>=0) {lo=a.lo*b.hi;hi=a.hi*b.lo;}
            else if(b.hi<=0) {lo=a.hi*b.hi;hi=a.lo*b.lo;}
            else {lo=a.lo*b.hi;hi=a.lo*b.lo;}
        } else if(b.lo>=0) {lo=a.lo*b.hi;hi=a.hi*b.hi;}
        else if(b.hi<=0) {lo=a.hi*b.lo;hi=a.lo*b.lo;}
        else {
            count(2);lo=std::min(a.lo*b.hi,a.hi*b.lo);
            hi=std::max(a.lo*b.lo,a.hi*b.hi);
        }
        count(2);return enclosure(lo,hi);
    }
    Interval square(Interval a) {
        if(zero(a)) return {0,0};
        count(2);double lo,hi;
        if(a.lo>=0) {lo=a.lo*a.lo;hi=a.hi*a.hi;}
        else if(a.hi<=0) {lo=a.hi*a.hi;hi=a.lo*a.lo;}
        else return {0,enclosure(0,std::max(a.lo*a.lo,a.hi*a.hi)).hi};
        auto r=enclosure(lo,hi);r.lo=std::max(0.,r.lo);return r;
    }
    Interval div(Interval a,Interval b) {
        if(b.lo<=0 && b.hi>=0) throw Failure(4);
        if(zero(a)) return {0,0};
        if(one(b)) return a;
        count(4);
        std::array<double,4> p={a.lo/b.lo,a.lo/b.hi,a.hi/b.lo,a.hi/b.hi};
        auto mm=std::minmax_element(p.begin(),p.end());return enclosure(*mm.first,*mm.second);
    }
};
struct Matrix {
    std::vector<std::map<I64,Interval>> rows;
    I64 budget, count=0, peak=0, created=0;
    Matrix(I64 n,I64 b): rows(n), budget(b) {if(n>b)throw Failure(1);}
    Interval get(I64 i,I64 j) const {
        auto it=rows[i].find(j);return it==rows[i].end()?Interval{0,0}:it->second;
    }
    void set(I64 i,I64 j,Interval v) {
        bool old=rows[i].count(j);
        if(Arithmetic::zero(v)) {
            if(old) {rows[i].erase(j);if(i!=j)rows[j].erase(i);--count;}
            return;
        }
        if(!std::isfinite(v.lo)||!std::isfinite(v.hi)||v.lo>v.hi)throw Failure(3,i,v);
        if(!old) {
            if(count>=budget)throw Failure(1,i);
            ++count;++created;peak=std::max(peak,count);
        }
        rows[i][j]=v;if(i!=j)rows[j][i]=v;
    }
    void remove(I64 i) {while(!rows[i].empty())set(i,rows[i].begin()->first,{0,0});}
};
std::array<I64,2> signature(Interval a,Interval b,Interval c,Arithmetic& ar,Interval& det) {
    det=ar.sub(ar.mul(a,c),ar.square(b));
    if(det.hi<0)return {1,1};
    if(det.lo>0) {
        auto tr=ar.add(a,c);
        if(tr.lo>0)return {2,0};
        if(tr.hi<0)return {0,2};
    }
    return {-1,-1};
}
void eliminate(Matrix& m,Arithmetic& ar,const std::vector<I64>& order,
               bool positive_only,Pivot* pivots,Summary& result) {
    I64 n=order.size();std::vector<I64> rank(n);std::vector<bool> active(n,true);
    for(I64 i=0;i<n;++i)rank[order[i]]=i;
    result={0,0,0,0,0,0,0,m.count};
    for(I64 k:order) {
        if(!active[k])continue;
        auto a=m.get(k,k);std::vector<I64> neighbors;
        for(auto const& pair:m.rows[k])if(pair.first!=k)neighbors.push_back(pair.first);
        auto sort_neighbors=[&](){std::sort(neighbors.begin(),neighbors.end(),[&](I64 i,I64 j){return rank[i]<rank[j];});};
        sort_neighbors();
        I64 l=-1;Interval b={0,0},c={0,0},det={0,0},inv={0,0};
        std::array<I64,2> sig;
        std::vector<std::array<Interval,2>> v,u;
        if(a.lo>0||a.hi<0) {
            sig=a.lo>0?std::array<I64,2>{1,0}:std::array<I64,2>{0,1};
            if(positive_only&&sig[1])throw Failure(5,k,a);
            inv=ar.div({1,1},a);
            for(I64 i:neighbors) {
                auto x=m.get(k,i);v.push_back({x,{0,0}});u.push_back({ar.mul(x,inv),{0,0}});
            }
        } else {
            for(I64 candidate:neighbors) {
                b=m.get(k,candidate);c=m.get(candidate,candidate);sig=signature(a,b,c,ar,det);
                if(sig[0]>=0){l=candidate;break;}
            }
            if(l<0)throw Failure(6,k,a);
            if(positive_only&&sig[1])throw Failure(5,k,a);
            for(auto const& pair:m.rows[l])if(pair.first!=k&&pair.first!=l)neighbors.push_back(pair.first);
            neighbors.erase(std::remove(neighbors.begin(),neighbors.end(),l),neighbors.end());
            std::sort(neighbors.begin(),neighbors.end());neighbors.erase(std::unique(neighbors.begin(),neighbors.end()),neighbors.end());sort_neighbors();
            auto invdet=ar.div({1,1},det),p=ar.mul(c,invdet),q=ar.mul(Arithmetic::neg(b),invdet),t=ar.mul(a,invdet);
            for(I64 i:neighbors) {
                auto x=m.get(k,i),y=m.get(l,i);v.push_back({x,y});
                u.push_back({ar.add(ar.mul(p,x),ar.mul(q,y)),ar.add(ar.mul(q,x),ar.mul(t,y))});
            }
        }
        result.width=std::max(result.width,static_cast<I64>(neighbors.size()));
        for(std::size_t ii=0;ii<neighbors.size();++ii)for(std::size_t jj=ii;jj<neighbors.size();++jj) {
            Interval correction;
            if(l<0)correction=ii==jj?ar.mul(ar.square(v[ii][0]),inv):ar.mul(u[ii][0],v[jj][0]);
            else correction=ar.add(ar.mul(u[ii][0],v[jj][0]),ar.mul(u[ii][1],v[jj][1]));
            auto i=neighbors[ii],j=neighbors[jj];m.set(i,j,ar.sub(m.get(i,j),correction));
        }
        m.remove(k);active[k]=false;if(l>=0){m.remove(l);active[l]=false;}
        pivots[result.count++]={k,l,sig[0],sig[1],a.lo,a.hi,b.lo,b.hi,c.lo,c.hi,det.lo,det.hi};
        result.positive+=sig[0];result.negative+=sig[1];
    }
    result.peak=m.peak;result.created=m.created;
}

// Appel privé : la couche Python valide tailles, CSR canonique, symétrie,
// finitude et permutation avant de fournir ces tableaux contigus.
extern "C" int vinkulum_inertia_binary64(
    I64 n,I64 s,I64 nr,const I64* dp,const I64* dj,const double* dv,
    const I64* mp,const I64* mj,const double* mv,const double* b,double gamma,
    const I64* permutation,I64 budget_ops,I64 budget_coeff,bool diagonal_mass,
    Pivot* mass_pivots,Pivot* kkt_pivots,Summary* summaries,I64* diagnostic,
    double* bad_interval,double* phases) noexcept {
    Arithmetic ar(budget_ops);I64 phase=0;
    auto start=std::chrono::steady_clock::now();
    auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    try {
        if(!environment_valid())throw Failure(8);
        std::vector<I64> order(permutation,permutation+n);
        if(diagonal_mass)summaries[0]={n,0,0,0,0,0,0,n};
        else {
            Matrix mass(n,budget_coeff);
            for(I64 i=0;i<n;++i)for(I64 k=mp[i];k<mp[i+1];++k)if(mj[k]>=i)mass.set(i,mj[k],Arithmetic::point(mv[k]));
            eliminate(mass,ar,order,true,mass_pivots,summaries[0]);
        }
        phases[0]=elapsed();phase=1;Matrix a(n+s,budget_coeff);
        for(I64 row=0;row<nr;++row)for(I64 ii=dp[row];ii<dp[row+1];++ii)for(I64 jj=ii;jj<dp[row+1];++jj) {
            auto i=dj[ii],j=dj[jj];a.set(i,j,ar.add(a.get(i,j),ar.mul(Arithmetic::point(dv[ii]),Arithmetic::point(dv[jj]))));
        }
        for(I64 i=0;i<n;++i) {
            for(I64 k=mp[i];k<mp[i+1];++k)if(mj[k]>=i)a.set(i,mj[k],ar.sub(a.get(i,mj[k]),ar.mul(Arithmetic::point(gamma),Arithmetic::point(mv[k]))));
            for(I64 j=0;j<s;++j)a.set(i,n+j,Arithmetic::point(b[i*s+j]));
        }
        phases[1]=elapsed()-phases[0];phase=2;
        for(I64 j=0;j<s;++j)order.push_back(n+j);
        eliminate(a,ar,order,false,kkt_pivots,summaries[1]);
        if(summaries[1].positive!=n||summaries[1].negative!=s)throw Failure(7);
        phases[2]=elapsed()-phases[0]-phases[1];diagnostic[0]=ar.operations;return 0;
    } catch(const Failure& e) {
        diagnostic[0]=ar.operations;diagnostic[1]=phase;diagnostic[2]=e.index;
        bad_interval[0]=e.interval.lo;bad_interval[1]=e.interval.hi;return static_cast<int>(e.code);
    } catch(...) {diagnostic[0]=ar.operations;diagnostic[1]=phase;return 9;}
}

// Entrée privée de contre-épreuve : auditer aussi l'arithmétique compilée.
extern "C" int vinkulum_interval_binary64(int operation,double alo,double ahi,
    double blo,double bhi,double* output) noexcept {
    try {
        if(!std::isfinite(alo)||!std::isfinite(ahi)||!std::isfinite(blo)||!std::isfinite(bhi)||alo>ahi||blo>bhi)return 10;
        if(!environment_valid())return 8;
        Arithmetic ar(100);Interval a={alo,ahi},b={blo,bhi},r;
        switch(operation) {
            case 0:r=ar.add(a,b);break;
            case 1:r=ar.sub(a,b);break;
            case 2:r=ar.mul(a,b);break;
            case 3:r=ar.div(a,b);break;
            case 4:r=ar.square(a);break;
            default:return 10;
        }
        output[0]=r.lo;output[1]=r.hi;return 0;
    } catch(const Failure& e) {return static_cast<int>(e.code);}
    catch(...) {return 9;}
}
