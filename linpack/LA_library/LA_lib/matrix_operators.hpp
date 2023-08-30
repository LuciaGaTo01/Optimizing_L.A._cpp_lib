//==================================================================================================
//  File:       matrix_operators.hpp
//
//  Summary:    This header defines the functions to perform different arithmetic operations
//              with matrices.
//==================================================================================================

#ifndef LA_LIB_MATRIX_OPERATORS_HPP
#define LA_LIB_MATRIX_OPERATORS_HPP

namespace LA_lib{

    template<typename T>
    class matrix;

    template <class M1, class M2, class M3, typename T>
    void add_matrix (M1 &m, M2 const &lhs, M3 const &rhs)
    {
        assert((lhs.rows() == rhs.rows() && lhs.columns() == rhs.columns()
                && "Both matrices must have the same size"));
                
        int rows = lhs.rows();
        int cols = lhs.columns();
        
        if constexpr (std::is_same<M1,matrix<T>>::value){
            m.resize(rows, cols);
        }

        for (int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j){
                m(i,j) = lhs(i,j) + rhs(i,j);
            }
        }
    }

    template <class M1, class M2, class M3, typename T>
    void sub_matrix (M1 &m, M2 const &lhs, M3 const &rhs)
    {
    	assert((lhs.rows() == rhs.rows() && lhs.columns() == rhs.columns()
                && "Both matrices must have the same size"));
                
        int rows = lhs.rows();
        int cols = lhs.columns();

        if constexpr (std::is_same<M1,matrix<T>>::value){
            m.resize(rows, cols);
        }

        for (int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j){
                m(i,j) = lhs(i,j) - rhs(i,j);
            }
        }
    }

    template <class M1, class M2, class M3, typename T>
    void mult_matrix (M1 &m, M2 const &lhs, M3 const &rhs)
    {
        assert((lhs.columns() == rhs.rows()
                && "First matrix columns must be equal to second matrix rows"));

        int rows = lhs.rows();
        int cols = rhs.columns();

        if constexpr (std::is_same<M1,matrix<T>>::value){
            m.resize(rows, cols);
        }

        T aux;
        for (int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j){
                aux = 0;
                for (int k = 0; k < rhs.rows(); ++k){
                    aux = aux + lhs(i,k) * rhs(k,j);
                }
                m(i,j) = aux;
            }
        }
    }


    template <class M1, class M2, typename T>
    void mult_matrix_times_scalar (M1 &m, M2 const &lhs, T const &rhs)
    {
        int rows = lhs.rows();
        int cols = lhs.columns();

        if constexpr (std::is_same<M1,matrix<T>>::value){
            m.resize(rows, cols);
        }

        for (int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j){
                m(i,j) = lhs(i,j) * rhs;
            }
        }
    }


    template <class M1, class M2, typename T>
    void divide_matrix_by_scalar (M1 &m, M2 const &lhs, T const &rhs)
    {
        assert (rhs != 0 && "Division by 0 is not defined");

        int rows = lhs.rows();
        int cols = lhs.columns();

        if constexpr (std::is_same<M1,matrix<T>>::value){
            m.resize(rows, cols);
        }

        for (int i = 0; i < rows; ++i){
            for (int j = 0; j < cols; ++j){
                m(i,j) = lhs(i,j)/rhs;
            }
        }
    }

    template <class M1, class M2>
    void plus_equal_matrix (M1 &lhs, M2 const &rhs)
    {
    	assert((lhs.rows() == rhs.rows() && lhs.columns() == rhs.columns()
                && "Both matrices must have the same size"));
        
        for (int i = 0; i < lhs.rows(); ++i){
            for (int j = 0; j < lhs.columns(); ++j){
                lhs(i,j) += rhs(i,j);
            }
        }
    }

    template <class M1, class M2>
    void minus_equal_matrix (M1 &lhs, M2 const &rhs)
    {
    	assert((lhs.rows() == rhs.rows() && lhs.columns() == rhs.columns()
                && "Both matrices must have the same size"));
        
        for (int i = 0; i < lhs.rows(); ++i){
            for (int j = 0; j < lhs.columns(); ++j){
                lhs(i,j) -= rhs(i,j);
            }
        }
    }

    template <class M1, typename T>
    void matrix_times_equal_scalar (M1 &lhs, T const &rhs)
    {
        for (int i = 0; i < lhs.rows(); ++i){
            for (int j = 0; j < lhs.columns(); ++j){
                lhs(i,j) *= rhs;
            }
        }
    }

    template <class M1, typename T>
    void matrix_divide_equal_scalar (M1 &lhs, T const &rhs)
    {
        assert (rhs != 0 && "Division by 0 is not defined");

        for (int i = 0; i < lhs.rows(); ++i){
            for (int j = 0; j < lhs.columns(); ++j){
                lhs(i,j) /= rhs;
            }
        }
    }

}


#endif //LA_LIB_MATRIX_OPERATORS_HPP
