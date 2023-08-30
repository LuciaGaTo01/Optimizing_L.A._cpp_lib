//==================================================================================================
//  File:       matrix_view.hpp
//
//  Summary:    This header defines the matrix_view class template, defined with the aim of
//              being able to access some submatrices inside a defined object.
//==================================================================================================

#ifndef LINEAR_ALGEBRA_MATRIX_VIEW_HPP
#define LINEAR_ALGEBRA_MATRIX_VIEW_HPP

namespace LA_lib {

    template <class M1, typename T1, class M2, typename T2>
    concept valid_matrix_view_op_type = (std::is_same<M1, matrix_view<T1>>::value
                                         && std::is_same<M2, matrix_view<T2>>::value)
                                        || (std::is_same<M1, matrix_transpose<matrix_view<T1>>>::value
                                            && std::is_same<M2, matrix_transpose<matrix_view<T2>>>::value)
                                        || (std::is_same<M1, matrix_view<T1>>::value
                                            && std::is_same<M2, matrix_transpose<matrix_view<T2>>>::value)
                                        || (std::is_same<M1, matrix_transpose<matrix_view<T1>>>::value
                                            && std::is_same<M2, matrix_view<T2>>::value);

    template <typename T>
    class matrix_view
    {
    private:
        matrix<T> &mview;
        int row_index;
        size_t nrows;
        int column_index;
        size_t ncols;

    public:

        using element_type = T;

// ==============================================================================
//                          MATRIX VIEW INITIALIZATION
// ==============================================================================

        template <typename N>
        matrix_view(matrix<T> &m, int ri, N rn, int ci, N cn):
                mview{m},
                row_index{ri}, nrows{static_cast<size_t>(rn)},
                column_index{ci}, ncols{static_cast<size_t>(cn)}
        {}

// ==============================================================================
//                            MATRIX VIEW DATA ACCESS
// ==============================================================================

        T operator()(int ri, int ci) const
        {
            return mview(row_index + ri, column_index + ci);
        }

        T &operator()(int ri, int ci)
        {
            return mview(row_index + ri, column_index + ci);
        }

        T operator()(int i) const
        {
            assert((nrows == 1 || ncols == 1 && "2D object cannot be accessed by single index"));

            if (nrows == 1) {
                return mview(1, column_index + i);
            } else {
                return mview(row_index + i, 1);
            }
        }

        T &operator()(int i)
        {
            assert((nrows == 1 || ncols == 1 && "2D object cannot be accessed by single index"));

            if (nrows == 1) {
                return mview(1, column_index + i);
            } else {
                return mview(row_index + i, 1);
            }
        }


// ==============================================================================
//                            MATRIX PROPERTIES
// ==============================================================================

        [[nodiscard]] constexpr bool
        is_empty() const noexcept {
            if (nrows == 0 && ncols == 0)
                return true;
            else
                return false;
        }

        [[nodiscard]] constexpr size_t
        rows() const noexcept { return nrows; }

        [[nodiscard]] constexpr size_t
        columns() const noexcept { return ncols; }

        [[nodiscard]] constexpr size_t
        size() const noexcept { return (nrows*ncols); }


// ==============================================================================
//                            MATRIX MODIFIERS
// ==============================================================================

        constexpr void
        swap_rows(int r1, int r2)
        {
            assert((r1 > 0 && r2 > 0 && "Invalid negative index"));

            assert((r1 < nrows && r2 < nrows && "Row index is out of range"));

            if (r1 != r2) {
                T aux;
                for (int i = 0; i < ncols; ++i) {
                    aux = this->operator()(r1, i);
                    this->operator()(r1, i) = this->operator()(r2, i);
                    this->operator()(r2, i) = aux;
                }
            }
        }

        constexpr void
        swap_columns(int c1, int c2)
        {
            assert((c1 > 0 && c2 > 0 && "Invalid negative index"));

            assert((c1 < nrows && c2 < nrows && "Column index is out of range"));

            if (c1 != c2) {
                T aux;
                for (int i = 0; i < nrows; ++i) {
                    aux = this->operator()(i, c1);
                    this->operator()(i, c1) = this->operator()(i, c2);
                    this->operator()(i, c2) = aux;
                }
            }
        }


// ==============================================================================
//                         ASSIGNMENT OPERATOR OVERLOADING
// ==============================================================================

        matrix_view<T> &operator=(matrix<T> const &rhs)
        {
            assert((nrows == rhs.rows() && ncols == rhs.columns()
                    && "Both matrices must have the same size"));
                    
            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = rhs(i,j);
                }
            }

            return *this;
        }

        matrix_view<T> &operator=(matrix_view<T> const &rhs)
        {
            assert((nrows == rhs.rows() && ncols == rhs.columns()
                    && "Both matrices must have the same size"));
                    
            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = rhs(i,j);
                }
            }

            return *this;
        }

        template <class M>
        matrix_view<T> &operator=(matrix_transpose<M> const &rhs)
        requires
        contain_same_element_type<typename M::element_type, T>
        {
            assert((nrows == rhs.rows() && ncols == rhs.columns()
                    && "Both matrices must have the same size"));
                    
            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = rhs(i,j);
                }
            }

            return *this;
        }


// ==============================================================================
//                      ARITHMETIC OPERATORS OVERLOADING
// ==============================================================================

        // --------------------------- ADDITION ---------------------------

        template <class M1, class M2>
        friend matrix<T> operator +(const M1 & lhs, const M2 & rhs)
        requires
        contain_same_element_type<typename M1::element_type, typename M2::element_type>
        and
        valid_matrix_view_op_type<M1, typename M1::element_type, M2, typename M2::element_type>
        {
            matrix<T> m;
            add_matrix<matrix<T>, M1, M2, T>(m, lhs, rhs);
            return m;
        }

        template <class M>
        matrix_view<T> & operator +=(const M &rhs)
        requires
        contain_same_element_type<T, typename M::element_type>
        {
            plus_equal_matrix<matrix_view<T>, M>(*this, rhs);
            return *this;
        }


        // --------------------------- SUBTRACTION ---------------------------

        template <class M1, class M2>
        friend matrix<T> operator -(const M1 & lhs, const M2 & rhs)
        requires
        contain_same_element_type<typename M1::element_type, typename M2::element_type>
        and
        valid_matrix_view_op_type<M1, typename M1::element_type, M2, typename M2::element_type>
        {
            matrix<T> m;
            sub_matrix<matrix<T>, M1, M2, T>(m, lhs, rhs);
            return m;
        }

        template <class M>
        matrix_view<T> & operator -=(const M &rhs)
        requires
        contain_same_element_type<T, typename M::element_type>
        {
            minus_equal_matrix<matrix_view<T>, M, T>(*this, rhs);
            return *this;
        }


        // --------------------------- MATRIX PRODUCT ---------------------------

        template <class M1, class M2>
        friend matrix<T> operator *(const M1 & lhs, const M2 & rhs)
        requires
        contain_same_element_type<typename M1::element_type, typename M2::element_type>
        and
        valid_matrix_view_op_type<M1, typename M1::element_type, M2, typename M2::element_type>
        {
            matrix<T> m;
            mult_matrix<matrix<T>, M1, M2, T>(m, lhs, rhs);
            return m;
        }

        template <class M>
        matrix_view<T> operator *=(const M & rhs)
        requires
        contain_same_element_type<typename M::element_type, T>
        {
            matrix<T> m;
            mult_matrix<matrix<T>, matrix_view<T>, M, T>(m, *this, rhs);
            *this = m;
            return *this;
        }


        // --------------------------- SCALAR PRODUCT ---------------------------

        template <class M>
        friend matrix<T> operator*(M const &lhs, T const &rhs)
        requires
        contain_same_element_type<typename M::element_type, T>
        and
        (std::is_same<M,matrix_view<T>>::value or std::is_same<M,matrix_transpose<matrix_view<T>>>::value)
        {
            matrix<T> m;
            mult_matrix_times_scalar<matrix<T>, M, T>(m, lhs, rhs);
            return m;
        }

        template <class M>
        friend matrix<T> operator*(T const &rhs, M const &lhs)
        requires
        contain_same_element_type<typename M::element_type, T>
        and
        (std::is_same<M,matrix_view<T>>::value or std::is_same<M,matrix_transpose<matrix_view<T>>>::value)
        {
            matrix<T> m;
            mult_matrix_times_scalar<matrix<T>, M, T>(m, lhs, rhs);
            return m;
        }

        matrix_view<T> & operator *=(const T &rhs)
        {
            matrix_times_equal_scalar<matrix_view<T>,T>(*this, rhs);
            return *this;
        }


        // --------------------------- DIVISION ---------------------------

        template <class M>
        friend matrix<T> operator/(M const &lhs, T const &rhs)
        requires
        contain_same_element_type<typename M::element_type, T>
        and
        (std::is_same<M,matrix_view<T>>::value or std::is_same<M,matrix_transpose<matrix_view<T>>>::value)
        {
            matrix<T> m;
            divide_matrix_by_scalar<matrix<T>, M, T>(m, lhs, rhs);
            return m;
        }

        matrix_view<T> & operator /=(const T &rhs)
        {
            matrix_divide_equal_scalar<matrix_view<T>,T>(*this, rhs);
            return *this;
        }


        // --------------------------- COMPARISON ---------------------------

        matrix_view<T> operator==(matrix_view<T> const &rhs)
        {
            if (nrows != rhs.rows() || ncols != rhs.columns())
                return false;

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    if (this->operator()(i,j) != rhs(i,j))
                        return false;
                }
            }

            return true;
        }

        matrix_view<T> operator!=(matrix_view<T> const &rhs)
        {
            if (nrows != rhs.rows() || ncols != rhs.columns())
                return true;

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    if (this->operator()(i,j) != rhs(i,j))
                        return true;
                }
            }

            return false;
        }


// ==============================================================================
//                               SUBELEMENT ACCESS
// ==============================================================================

        matrix_transpose<matrix_view<T>> t(){
            return matrix_transpose<matrix_view<T>>{*this};
        }

    };
}

#endif //LINEAR_ALGEBRA_MATRIX_VIEW_HPP
