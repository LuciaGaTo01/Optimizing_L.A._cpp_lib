//==================================================================================================
//  File:       matrix_view.hpp
//
//  Summary:    This header defines the matrix_view class template, defined with the aim of
//              being able to access some submatrices inside a defined object.
//==================================================================================================

#ifndef LINEAR_ALGEBRA_MATRIX_VIEW_HPP
#define LINEAR_ALGEBRA_MATRIX_VIEW_HPP

namespace LA_lib {

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
            return mview.data()[((row_index + ri) * mview.columns()) + column_index + ci];
        }

        T &operator()(int ri, int ci)
        {
            return mview.data()[((row_index + ri) * mview.columns()) + column_index + ci];
        }

        T operator()(int i) const
        {
            assert(nrows == 1 || ncols == 1 && "2D object cannot be accessed by single index");

            if (nrows == 1) {
                return mview.data()[column_index + i];
            } else {
                return mview.data()[row_index + i];
            }
        }

        T &operator()(int i)
        {
            assert(nrows == 1 || ncols == 1 && "2D object cannot be accessed by single index");

            if (nrows == 1) {
                return mview.data()[column_index + i];
            } else {
                return mview.data()[row_index + i];
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
            assert(r1 > 0 && r2 > 0 && "Invalid negative index");

            assert(r1 < nrows && r2 < nrows && "Row index is out of range");

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
            assert(c1 > 0 && c2 > 0 && "Invalid negative index");

            assert(c1 < nrows && c2 < nrows && "Column index is out of range");

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

        // --------------------------- GENERAL ASSIGNMENT ---------------------------

        template<class M>
        matrix_view<T> &operator=(M const &rhs)
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


        // --------------------------- ADDITION ---------------------------

        template <class M>
        matrix_view<T> & operator +=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            plus_equal_matrix<matrix_view<T>, M>(*this, rhs);
            return *this;
        }


        // --------------------------- SUBTRACTION ---------------------------

        template <class M>
        matrix_view<T> & operator -=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            minus_equal_matrix<matrix_view<T>, M>(*this, rhs);
            return *this;
        }

        matrix<T> operator-()
        {
            matrix<T> m{nrows, ncols};
            for (int i = 0; i < nrows; ++i) {
                for (int j = 0; j < ncols; ++j) {
                    m(i, j) = -this->operator()(i, j);
                }
            }
            return m;
        }


        // --------------------------- MATRIX PRODUCT ---------------------------

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

        matrix_view<T> & operator *=(const T &rhs)
        {
            matrix_times_equal_scalar<matrix_view<T>,T>(*this, rhs);
            return *this;
        }


        // --------------------------- DIVISION ---------------------------

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
