/// \author Stefano Veroni
/// \file VLabs

#ifndef MYCFUNCTIONS_H
#define MYCFUNCTIONS_H


#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>


using namespace std;


//_____________________________________________________________________
//
//---------------------------------------------------------------------
// VECTORS ADDITIONAL FUNCTIONS
//---------------------------------------------------------------------
//_____________________________________________________________________
//From https://thispointer.com/cpp-vector-print-all-elements/
template <typename T>
inline
void print_vector(const std::vector<T> & vec, 
                  std::string sep=" , ",
                  bool brackets = true)
    /** @brief rints given vector, with given separator.
    *
    *  @param vec The vector to print. Type: std::vector<whatever>.
    *  @param sep Separator to use between entries. Type: Char. Default " , ".
    *  @param brackets Start and end with brackets?. Type: bool. Default true.
    * 
    *  @return Void.
    **/
    {   
        int size = vec.size();
        if (brackets) {std::cout<<"\n{ ";}
        for(int index = 0; index < size ; index++)
        {   
            if (index == 0) {std::cout<<vec[index];}
            else {std::cout<<sep<<vec[index];}
        }
        if (brackets) {std::cout<<" }";}
        std::cout<<std::endl;  
}


template <typename T>
inline
void print_vector(std::vector<vector<T>> vec, 
                  std::string sep=" , ",
                  bool brackets = true)
    /** @brief rints given vector, with given separator.
    *
    *  @param vec The vector to print. Type: std::vector<vector<whatever>>.
    *  @param sep Separator to use between entries. Type: Char. Default " , ".
    *  @param brackets Start and end with brackets?. Type: bool. Default true.
    * 
    *  @return Void.
    **/
    {   int size = vec.size();
        if (brackets) {std::cout<<"\n{ ";}
        for(int index = 0; index < size ; index++)
        {
            std::vector<T> vind = vec[index];
            std::cout<<"{";
            for (int j = 0; j < vind.size(); j++)
            {
                if (j == 0) {std::cout<<vind[0];}
                else {std::cout<<sep<<vind[j];}
            }
            std::cout<<"}";

            if (index != size -1) {cout<<sep;}
        }
        if (brackets) {std::cout<<" }";}
        std::cout<<std::endl;  
}


template <typename T>
inline
T myvecsum(vector <T> v)
    {return accumulate(v.begin(),v.end(), .0);}

template <typename T, typename F>
inline
double myvecsum(vector <T> v, F func)
{
    double sum = 0;
    for (int i = 0; i < v.size(); i++) sum += func(v[i]);
    return sum;
}


template <typename T1>
inline
T1 vecmax(vector<T1> v)
    {return *(max_element(v.begin(), v.end()));}

template <typename T1>
inline
T1 vecmin(vector<T1> v)
    {return *(min_element(v.begin(), v.end()));}


template <typename T1, typename T2>
int argmax(std::vector<T1, T2> const& v, int begin = 0, int end = 0) 
/** @brief Gets index of maximum in vector in interval [begin, end)
 *
 *  @param v     Vector where to look in.
 *  @param begin Starting point of range. Type int. 
 *               Default 0.
 *  @param end   Endng point of range. Type int.
 *               Default 0 (meaning v.end() is end of range). 
 * 
 *  @return Index of maximum.
 **/
{
    auto b = v.begin() + begin;
    auto e = (end == 0) ? v.end() : v.begin()+end;
    return static_cast<int>(std::distance(v.begin(), 
                                          max_element(b, e)));
}

template <typename T1, typename T2>
int argmin(std::vector<T1, T2> const& v, int begin = 0, int end = 0) 
{
/** @brief Gets index of minimum in vector in interval [begin, end)
 *
 *  @param v     Vector where to look in.
 *  @param begin Starting point of range. Type int. 
 *               Default 0.
 *  @param end   Endng point of range. Type int.
 *               Default 0 (meaning v.end() is end of range). 
 * 
 *  @return Index of minimum.
 **/
    auto b = v.begin() + begin;
    auto e = (end == 0) ? v.end() : v.begin()+end;
    return static_cast<int>(std::distance(v.begin(), 
                                          min_element(b, e)));
}


template <typename T>
inline
vector <int> getIndexes(vector<T> v, T x)
/** @brief Gets indexes where given vector's entries have value i.
 *
 *  @param v Vector to look in. Type: std::vector<whatever>.
 *  @param x Value to look for. Type: whatever. 
 * 
 *  @return Vector of indexes. Type: vector <int>
 **/
{
    vector <int> sol;
    auto it = v.begin();
    while (it != v.end()) 
    {
        auto it2 = find(it+1, v.end(), x);
        if (it2 != v.end()) 
        {
        int index = it2 - v.begin();
        sol.push_back(index);
        }
        it = it2;
    }
    return sol;
}

template <typename T1, typename T2>
inline
bool contains(vector <T1> v, T2 x)
    {return std::find(v.begin(), v.end(), x) != v.end();}


template <typename T1, typename T2>
inline
vector <T1> getAwhereB(vector<T1> v1, vector<T2> v2, T2 x)
/** @brief Gets values of vector v1 at indexes where vector v2 = x.
 *
 *  @param v1 Vector to take values from. Type: std::vector<whatever>.
 *  @param v2 Vector where to look for. Type: std::vector<whatever2>. 
 *  @param x  Value to look for. Type: whatever2. 
 * 
 *  @return Values of vector v1 at indexes where vector v2 = x.
 **/
{   
    if (v1.size() != v2.size())
    {
        throw std::invalid_argument( "Vectors Have Different Sizes" );
    } 
    vector <T1> sol;
    auto begin = v2.begin();
    auto end = v2.end();
    auto it = begin-1;
    while (it != end) 
    {
        auto it2 = find(it+1, end, x);
        if (it2 != end) 
        {
        int index = it2 - begin;
        sol.push_back(v1[index]);
        }
        it = it2;
    }
    return sol;
}

template <typename T1, typename T2>
inline
vector<vector<T1>> getAwhereB(vector<vector<T1>> V, vector<T2> u, T2 x)
/** @brief Gets values of vectors v0=V[0], v1=V[1], etc... at indexes 
 *         where vector U = x.
 *
 *  @param V  Vector of many vectors to take values from. 
 *            Type: std::vector< vector<whatever> >.
 *  @param U Vector where to look for. Type: std::vector<whatever2>. 
 *  @param x  Value to look for. Type: whatever2. 
 * 
 *  @return Vector of vectors with values corrsponding to indexes where
 *          vector U = x.
 **/
{   
    for (int i = 0; i < V.size(); i++)
    {
        if (V[i].size() != u.size())
        {
            throw std::invalid_argument
            ( "Vector " + to_string(i) + " in V and Vector U have Different Sizes" );
        } 
    }
    vector< vector<T1> > sol (V.size()); //one vector per each vector in V
    auto begin = u.begin();
    auto end = u.end();
    auto it = begin-1;
    while (it != end) 
    {
        auto it2 = find(it+1, end, x);
        if (it2 != end) 
        {
        int index = it2 - begin;
        for (int i = 0; i < V.size(); i++)
            {
                sol[i].push_back(V[i][index]);
            }
        }
        it = it2;
    }
    return sol;
}

template <typename T1, typename T2, typename F>
vector <T1> getAwhereB(vector<T1> v1, vector<T2> v2, F f)
/** @brief Gets values of vector v1 at indexes where v2 satisfies f.
 *
 *  @param v1 Vector to take values from. Type: std::vector<whatever>.
 *  @param v2 Vector where to look for. Type: std::vector<whatever2>. 
 *  @param f  Function taking as argument whatever2 type and 
 *            returning bool. 
 * 
 *  @return Values of vector v1 at indexes where vector v2 satisfies f.
 **/
{   
    if (v1.size() != v2.size())
    {
        throw std::invalid_argument( "Vectors Have Different Sizes" );
    } 
    vector <T1> sol;
    auto begin = v2.begin();
    auto end = v2.end();
    auto it = begin-1;
    for (int ind = 0; ind < v2.size(); ind++) 
        {if ( f(v2[ind]) == 1 ) sol.push_back(v1[ind]);}
    return sol;
}

template <typename T1, typename T2, typename F>
vector <vector<T1>> getAwhereB(vector<vector<T1>> v1, vector<T2> v2, F f)
/** @brief Gets values of vector v1 at indexes where v2 satisfies f.
 *
 *  @param v1 Vector to take values from. Type: std::vector<whatever>.
 *  @param v2 Vector where to look for. Type: std::vector<whatever2>. 
 *  @param f  Function taking as argument whatever2 type and 
 *            returning bool. 
 * 
 *  @return Values of vector v1 at indexes where vector v2 satisfies f.
 **/
{   
    if (v1.size() != v2.size())
    {
        throw std::invalid_argument( "Vectors Have Different Sizes" );
    } 
    vector <T1> sol;
    auto begin = v2.begin();
    auto end = v2.end();
    auto it = begin-1;
    for (int ind = 0; ind < v2.size(); ind++) 
        {if ( f(v2[ind]) == 1 ) sol.push_back(v1[ind]);}
    return sol;
}


template <typename T1, typename T2>
vector <T1> sortAwithB(vector <T1> A, vector <T2> B, bool reverse = false)
{
    struct mypair
    {
        T1 a;
        T2 b;
    };

    vector <mypair> zipped (A.size());
    for (int i = 0; i<A.size(); i++)
    {
        zipped[i] = {A[i], B[i]};
    }

    if (not reverse)
    sort(zipped.begin(), zipped.end(), 
         [](const auto& i, const auto& j){ return i.b < j.b; } );
    else
    sort(zipped.begin(), zipped.end(), 
         [](const auto& i, const auto& j){ return i.b > j.b; } );
    
    vector <T1> sortA (A.size());
    for (int i = 0; i<A.size(); i++)
    {
        sortA[i] = zipped[i].a;
    }
    
    return sortA;
}


template <typename T1>
vector <T1> sum_portions(vector <T1> v, int n)
/** @brief Sums first, second, ..., nth portion of vector 
 *         with each other.
 *
 *  @param v Vector. Type vector <whatever>.
 *  @param n Number of portions into which divide v, then sum
 *           on each other. Type int.

 *  @return Vector of v.size()/n elements. 
 **/
{
    if (v.size()%n) 
    { 
        throw std::invalid_argument("Vector's size not multiple of n");
    }
    int m = v.size()/n;
    vector <T1> sol (m); 
 
    for (int i = 0; i < v.size(); i++) 
    {
        int index = i%m;
        sol[index] += v[i];
    }
    return sol;
}


template <typename T1, typename T2>
inline
double minAbsDiff (vector<T1> v1, vector<T2> v2)
/** @brief Minimum absolute difference between elements of v1 and v2.
 *
 *  @param v1 vector 1 of numbers
 *  @param v2 vector 2 of numbers
 * 
 *  @return Minimum absolute difference between elements of v1 and v2.
 *          Type: double.
 **/
{
    double min = INFINITY;  
    for (int i=0; i<v1.size(); i++)
    {
        for (int j=0; j<v2.size(); j++)
        {
            double diff = abs( v1[i] - v2[j] );
            min = (diff < min) ? diff : min;
        }
    }
    return min;
}


//https://gist.github.com/lorenzoriano/5414671
template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

//_____________________________________________________________________
//
//---------------------------------------------------------------------
// PROPER VECTOR MATH FUNCTIONS
//---------------------------------------------------------------------
//_____________________________________________________________________

template <typename T1>
inline
double mymean(vector <T1> x)
    {return (double)myvecsum(x) / x.size();}

template <typename T1, typename F>
inline
double mymean(vector <T1> x, F func)
    {return (double)myvecsum(x, func) / x.size();}

template <typename T1, typename T2>
inline
double mymean(vector <T1> x, vector <T2> w = {1})
{   
    if (x.size() != w.size())
    {
        throw std::invalid_argument
        ("vector and weight have different sizes: "
        +to_string(x.size()) +" and " 
        +to_string(w.size()) + "!");
    }
    T2 sumw = myvecsum(w);
    double mean = 0;
    for (int i = 0; i < x.size(); i++)
        {mean += (double) x[i] * w[i] / sumw;}
    return mean;
}

template <typename T1, typename T2, typename F>
inline
double mymean(vector <T1> x, F func, vector <T2> w = {1})
{   
    if (w.size () == 1) return mymean(x, func);

    else if (x.size() != w.size())
    {
        throw std::invalid_argument
        ("vector and weight have different sizes: "
        +to_string(x.size()) +" and " 
        +to_string(w.size()) + "!");
    }

    T2 sumw = myvecsum(w);
    double mean = 0;
    for (int i = 0; i < x.size(); i++)
        {mean += (double) func(x[i]) * w[i];}
    return (double) mean / sumw;
}


template <typename T>
inline
T cumulative (vector <T> v)
{
    vector <double> cum (N);
    T tot = myvecsum(v);
    v[0] /= tot;
    partial_sum(v.begin(), v.end(), cum.begin(), 
               [](double a, double b)
                {return (double)a + b;});
    return cum;
}

template <typename T>
inline
double cumulative (vector <T> v, bool norm = true)
{
    vector <double> cum (N);
    T tot = myvecsum(v);
    v[0] /= tot;
    if (norm) auto lambda = [tot](double a, double b)
                            {b/=tot; return (double)a + b;};
    else      auto lambda = [](double a, double b)
                            {return (double)a + b;};
    partial_sum(v.begin(), v.end(), cum.begin(), lambda);
    return cum;
}


template <typename T>
inline
int get_sign (T x)
{
    return (x > 0) - (x < 0);
}

struct Point
        {double x, y;};

template <typename T>
//https://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/ 
vector <vector<T>> ConvexHull(vector<T> xvec, vector<T> yvec)
{
    const int N = xvec.size();
    if (N <= 3)
    {
        vector <vector<T>> sol (N);
        for (int i = 0; i<N; i++)
            {sol[i] = {xvec[i], yvec[i]};}
        return sol;
    }
    else
    {
        vector <Point> points (N);
        for (int i = 0; i<N; i++)
            {points[i] = {xvec[i], yvec[i]};}
        return ConvexHull(points, 0);
    }
}

template <typename T>
//https://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/ 
vector <vector<T>> ConvexHull(vector <vector<T>> positions)
{
    const int N = positions.size();
    if (N <= 3)
        {return positions;}
    else
    {
        vector <Point> points (N);
        for (int i = 0; i<N; i++)
            {points[i] = {positions[i][0], positions[i][1]};}
        return ConvexHull(points, 0);
    }
}


template <typename T>
inline
vector <vector<T>> StarHull(vector<T> xvec, vector<T> yvec, 
                            T centrex, T centrey)
/** @brief Sort points in order about centre
 * 
 *  @return {sortedx, sortedy}
 **/
{
    vector <double> phi (xvec);
    for (int i = 0; i<xvec.size(); i++)
        {
            double dyi = yvec[i] - centrey;
            double dxi = xvec[i] - centrex;
            phi[i] = atan2(dyi, dxi);
        }
    vector <vector<T>> sorted = {sortAwithB(xvec, phi), 
                                 sortAwithB(yvec, phi)};
    return sorted;
}


template <typename T>
inline
double Shoelace (vector<T> xvec, vector<T> yvec)
/** @brief Get Area of figure given ordered set of x and y coords.
 *
 *  @param xvec 
 *  @param yvec
 * 
 *  @return Area
 **/
{
    // Area
    double A = 0.0;
 
    // Calculate value of shoelace formula
    int n = xvec.size();
    int j = n - 1;
    for (int i = 0; i < n; i++)
    {
        A += (xvec[j] + xvec[i]) * (yvec[j] - yvec[i]);
        j = i;  // j is previous vertex to i
    }
 
    // Return absolute value
    return abs(A / 2.0);
}

template <typename T>
inline
double Shoelace (vector<vector<T>> posvec)
/** @brief Get Area of figure given a vector of points in 2D.
 *
 *  @param posvec Vector of 2D vectors. 
 * 
 *  @return Area
 **/
{
    // Area
    double A = 0.0;
    int n = posvec.size();
    vector <double> xvec (n), yvec(n);

    for (int i = 0; i < n; i++)
    {
        xvec[i] = posvec[i][0];
        yvec[i] = posvec[i][1];
    }
    // Calculate value of shoelace formula
    int j = n - 1;
    for (int i = 0; i < n; i++)
    {
        A += (xvec[j] + xvec[i]) * (yvec[j] - yvec[i]);
        j = i;  // j is previous vertex to i
    }
 
    // Return absolute value
    return abs(A / 2.0);
}




/////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------
// ADDITIONAL STUFF YOU ARE LIKELY TO NEVER USE 
// BUT THEY ARE USED INSIDE SOME OF ABOVE FUNCTIONS
/////////////////////////////////////////////////////////////////////////

// A utility function to find next to top in a stack
inline
Point nextToTop(stack<Point> &S)
    {
        Point p = S.top();
        S.pop();
        Point res = S.top();
        S.push(p);
        return res;
    }

// A utility function to swap two points
inline
void swap(Point &p1, Point &p2)
    {
        Point temp = p1;
        p1 = p2;
        p2 = temp;
    }

// A utility function to return square of distance
// between p1 and p2
inline
double distSq(Point p1, Point p2)
    {
        return (p1.x - p2.x)*(p1.x - p2.x) +
            (p1.y - p2.y)*(p1.y - p2.y);
    }

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
inline
double orientation(Point p, Point q, Point r)
    {
        //(y2 - y1)*(x3 - x2) - (y3 - y2)*(x2 - x1)
        double val = (q.y - p.y) * (r.x - q.x) -
                  (q.x - p.x) * (r.y - q.y);
    
        if (val == 0) return 0;  // collinear
        return (val > 0)? 1: 2;  // clock or counterclock wise
    }

inline
double get_phi(Point p, Point p0)
{
    double dy = p.y - p0.y;
    double dx = p.x - p0.x;
    double phi = atan2(dy, dx);//+ M_PI;
    return phi;
}

template <typename T>
inline
double get_phi(vector<T> p, vector<T> p0)
{
    double dy = p[1] - p0[1];
    double dx = p[0] - p0[0];
    double phi = atan2(dy, dx);// + M_PI +1e-10;
    return phi;
}

//convex hull of a set of n points.
vector <vector<double>> ConvexHull(vector <Point> points, int n = 0)
    {
    // Find the bottommost point
    n = points.size();
    double ymin = points[0].y;
    int min = 0;
    for (int i = 1; i < n; i++)
    {
        double y = points[i].y;
        // Pick the bottom-most or chose the left
        // most point in case of tie
        if ((y<ymin) || (y==ymin && points[i].x < points[min].x))
        {   
            ymin = points[i].y;
            min = i;
        }
    }
    
    // Place the bottom-most point at first position
    Point p0 = points[min];
    points[min] = points[0];
    points[0] = p0;
    //cout<<points[min].x<<", "<<points[min].y<<endl;
    

    // Sort n-1 points with respect to the first point.
    // A point p1 comes before p2 in sorted output if p2
    // has larger polar angle (in counterclockwise
    // direction) than p1
    // If two or more points make same angle with p0,
    // Remove all but the one that is farthest from p0
    // Remember that, in above sorting, our criteria was
    // to keep the farthest point at the end when more than
    // one points have same angle.
    vector <double> phis (n);
    for (int i = 1; i<n; i++)
        {phis[i] = get_phi(points[i], p0);}
    points = sortAwithB(points, phis);
    sort(phis.begin(), phis.end());
    //for (int i = 0; i<n; i++)
    //{cout<<(points[i].y-points[0].y)/(points[i].x-points[0].x)<<endl;}
    //print_vector(phis);

    for (int i = 0; i<n-1;)
    {
        if (i+1>=phis.size()) {break;}
        //cout<<i<<": "<<phis[i]<<" vs "<<phis[i+1]<<endl;
        //print_vector(phis);
        if (phis[i] == phis[i+1])
        {
            phis.erase(phis.begin() + i);
            int j = (distSq(p0, points[i])>distSq(p0, points[i+1]))?
                    i+1 : i;
            /*cout<<endl;
            cout<<points[i].x<<", "<<points[i].y<<endl;
            cout<<points[i+1].x<<", "<<points[i+1].y<<endl;
            cout<<endl;*/
            points.erase(points.begin() + j);
        }
        else {i++;}
    }
    int m = points.size();
    //for (int i = 0; i<m; i++)cout<<points[i].x<<", "<<points[i].y<<endl;
    
    // Create an empty stack and push first three points
    // to it.
    stack<Point> S;
    S.push(points[0]);
    S.push(points[1]);
    S.push(points[2]);
    
    // Process remaining n-3 points
    for (int i = 3; i < m; i++)
    {
        // Keep removing top while the angle formed by
        // nextToTop(S), S.top(), points[i] makes
        // a non-left turn (non-counterclockwise path)
        while (S.size()>2 && 
                orientation(nextToTop(S), S.top(), points[i]) != 2)
            {
                /*
                cout<<orientation(nextToTop(S), S.top(), points[i])<<endl;
                cout<<nextToTop(S).x<<", "<<nextToTop(S).y<<endl;
                cout<<S.top().x<<", "<<S.top().y<<endl;
                cout<<points[i].x<<", "<<points[i].y<<endl;
                */
                S.pop();
            }
        S.push(points[i]);
    }
    
    // Now return solution as vector
    vector <vector<double>> sol;
    while (!S.empty())
    {
        Point p = S.top();
        sol.push_back({p.x, p.y});
        //cout << "(" << p.x << ", " << p.y <<")" << endl;
        S.pop();
    }
    reverse(sol.begin(), sol.end());
    return sol;
    }



#endif /* MYCFUNCTIONS_H */