#include <iostream>
#include <vector>
#include <initializer_list>
#include <typeinfo>
#include "linearAlgebra.hpp"


int main (){

    std::vector<int> v(2,3);
    LA_lib::matrix<int> m(2,1);

    m = {1,2};

    auto mm = m;

    m.resize_columns(2);

    // auto t = typeid(v + v);

    if (t==typeid(int)){
        std::cout << "Si soy \n\n";
    }
    // m(0,0) = 2;
    std::cout << m.columns() << "\n";
    std::cout << m(0,0) << "\n\n";

    std::cout << mm.rows() << "\n";
    std::cout << mm(0,1) << "\n";

    auto c = v;

    for (int i = 0; i < 2; ++i){
        c[i] = -v[i];
    }

    for (auto k: c){
        std::cout << k << " ";
    }
    std::cout << '\n';
}