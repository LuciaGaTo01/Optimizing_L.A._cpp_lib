//==================================================================================================
//  File:       checker.hpp
//
//  Summary:    This header define some functions that are used for checking some necessary
//              conditions for the library. It also provide PRINT() function for the matrices.
//==================================================================================================

#ifndef LA_LIB_CHECKER_HPP
#define LA_LIB_CHECKER_HPP

namespace LA_lib{

    template <typename T>
    const bool check_equal_size_rows (std::initializer_list<std::initializer_list<T>> const &init){
        T aux = init.begin()->size();

        for (auto v:init){
            if (v.size() != aux){
                return false;
            }
        }

        return true;

    }

    template <typename T>
    const bool is_list_empty (std::initializer_list<std::initializer_list<T>> const &init){

        if (std::empty(*init.begin()))
        {
            return true;
        }
        return false;

    }

    template <class T>
    void PRINT (T &mat){
        std::cout << "\n\nSize: " << mat.rows() << "x" << mat.columns() << "\n";
        std::cout << "Total number of elements: " << mat.size() << "\n";
        std::cout << "---------------------------------------------------------\n";
        for (int i = 0; i < mat.rows(); ++i){
            for (int j = 0; j < mat.columns(); ++j){
                std::cout << mat(i,j) << "\t";
            }
            std::cout << "\n";
        }
    }

    template <class T>
    void PRINT (T &&mat){
        std::cout << "\n\nSize: " << mat.rows() << "x" << mat.columns() << "\n";
        std::cout << "Total number of elements: " << mat.size() << "\n";
        std::cout << "---------------------------------------------------------\n";
        for (int i = 0; i < mat.rows(); ++i){
            for (int j = 0; j < mat.columns(); ++j){
                std::cout << mat(i,j) << "\t";
            }
            std::cout << "\n";
        }
    }
}



#endif //LA_LIB_CHECKER_HPP
