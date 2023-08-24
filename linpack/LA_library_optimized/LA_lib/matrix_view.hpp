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
            return mview.data()[((row_index + ri) * mview.columns()) + column_index + ci];
        }

        T &operator()(int ri, int ci)
        {
            return mview.data()[((row_index + ri) * mview.columns()) + column_index + ci];
        }

        T operator()(int i) const
        {
            if (nrows > 1 && ncols > 1)
                throw std::invalid_argument("2D object cannot be accessed by single index");

            if (nrows == 1) {
                return mview.data()[column_index + i];
            } else {
                return mview.data()[row_index + i];
            }
        }

        T &operator()(int i)
        {
            if (nrows > 1 && ncols > 1)
                throw std::invalid_argument("2D object cannot be accessed by single index");

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
            if (r1 < 0 || r2 < 0)
                throw std::invalid_argument("Invalid negative index");

            if (r1 >= nrows || r2 >= nrows)
                throw std::out_of_range("Row index is out of range");

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
            if (c1 < 0 || c2 < 0)
                throw std::invalid_argument("Invalid negative index");

            if (c1 >= ncols || c2 >= ncols)
                throw std::out_of_range("Column index is out of range");

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
            if (nrows != rhs.rows() || ncols != rhs.columns())
                throw std::domain_error("Both matrices must have the same size");

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = rhs(i,j);
                }
            }

            return *this;
        }

        matrix_view<T> &operator=(matrix_view<T> const &rhs)
        {
            if (nrows != rhs.rows() || ncols != rhs.columns())
                throw std::domain_error("Both matrices must have the same size");

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
            if (nrows != rhs.rows() || ncols != rhs.columns())
                throw std::domain_error("Both matrices must have the same size");

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
        matrix_view<T> & operator =(const matrix_template_addition<M1, M2> &mv)
        requires
                contain_same_element_type<T, typename M1::element_type>
        {
            add_matrix<matrix_view<T>, M1, M2, T>(*this, mv.lhs, mv.rhs);
            return *this;
        }

        template <class M>
        matrix_view<T> & operator +=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            plus_equal_matrix<matrix_view<T>, M>(*this, rhs);
            return *this;
        }

        // The following implementations allow to compute multiple operations in an efficient way.
        // No temporal objects are generated to obtain the result.

        template <class M1, class M2>
        matrix_view<T> & operator +=(const matrix_template_addition<M1, M2> &mv)
        requires
        contain_same_element_type<typename M1::element_type, typename M2::element_type>
        and
        contain_same_element_type<T, typename M1::element_type>
        {
            if (mv.lhs.rows() != mv.rhs.rows() || mv.lhs.columns() != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");

            if (nrows != mv.rhs.rows() || ncols != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");
                
            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) + (mv.lhs(i,j) + mv.rhs(i,j));
                }
            }

            return *this;
        }

        template <class M1, class M2>
        matrix_view<T> & operator +=(const matrix_template_subtraction<M1, M2> &mv)
        requires
                contain_same_element_type<typename M1::element_type, typename M2::element_type>
                and
                contain_same_element_type<T, typename M1::element_type>
        {
            if (mv.lhs.rows() != mv.rhs.rows() || mv.lhs.columns() != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");

            if (nrows != mv.rhs.rows() || ncols != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) + (mv.lhs(i,j) - mv.rhs(i,j));
                }
            }

            return *this;
        }

        template <class M>
        matrix_view<T> & operator +=(const matrix_template_scalar_product<M, T> &mv)
        requires
                contain_same_element_type<T, typename M::element_type>
        {

            if (nrows != mv.mat_view.rows() || ncols != mv.mat_view.columns())
                throw std::domain_error("All matrices must have the same size");

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) + (mv.mat_view(i,j) * mv.factor);
                }
            }

            return *this;
        }

        template <class M>
        matrix_view<T> & operator +=(const matrix_template_scalar_division<M, T> &mv)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            if (nrows != mv.mat_view.rows() || ncols != mv.mat_view.columns())
                throw std::domain_error("All matrices must have the same size");

            if (mv.factor == 0)
                throw std::overflow_error("Division by 0 is not defined");

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) + (mv.mat_view(i,j)/mv.factor);
                }
            }

            return *this;
        }


        // --------------------------- SUBTRACTION ---------------------------

        template <class M1, class M2>
        matrix_view<T> & operator =(const matrix_template_subtraction<M1, M2> &mv)
        requires
                contain_same_element_type<T, typename M1::element_type>
        {
            sub_matrix<matrix_view<T>,M1, M2, T>(*this, mv.lhs, mv.rhs);
            return *this;
        }

        template <class M>
        matrix_view<T> & operator -=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            minus_equal_matrix<matrix_view<T>, M>(*this, rhs);
            return *this;
        }

        // The following implementations allow to compute multiple operations in an efficient way.
        // No temporal objects are generated to obtain the result.

        template <class M1, class M2>
        matrix_view<T> & operator -=(const matrix_template_addition<M1, M2> &mv)
        requires
                contain_same_element_type<typename M1::element_type, typename M2::element_type>
                and
                contain_same_element_type<T, typename M1::element_type>
        {
            if (mv.lhs.rows() != mv.rhs.rows() || mv.lhs.columns() != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");

            if (nrows != mv.rhs.rows() || ncols != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) - (mv.lhs(i,j) + mv.rhs(i,j));
                }
            }

            return *this;
        }

        template <class M1, class M2>
        matrix_view<T> & operator -=(const matrix_template_subtraction<M1, M2> &mv)
        requires
        contain_same_element_type<typename M1::element_type, typename M2::element_type>
        and
        contain_same_element_type<T, typename M1::element_type>
        {
            if (mv.lhs.rows() != mv.rhs.rows() || mv.lhs.columns() != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");

            if (nrows != mv.rhs.rows() || ncols != mv.rhs.columns())
                throw std::domain_error("All matrices must have the same size");
                
            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) - (mv.lhs(i,j) - mv.rhs(i,j));
                }
            }

            return *this;
        }

        template <class M>
        matrix_view<T> & operator -=(const matrix_template_scalar_product<M, T> &mv)
        requires
                contain_same_element_type<T, typename M::element_type>
        {

            if (nrows != mv.mat_view.rows() || ncols != mv.mat_view.columns())
                throw std::domain_error("All matrices must have the same size");

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) - (mv.mat_view(i,j) * mv.factor);
                }
            }

            return *this;
        }

        template <class M>
        matrix_view<T> & operator -=(const matrix_template_scalar_division<M, T> &mv)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            if (nrows != mv.mat_view.rows() || ncols != mv.mat_view.columns())
                throw std::domain_error("All matrices must have the same size");

            if (mv.factor == 0)
                throw std::overflow_error("Division by 0 is not defined");

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = this->operator()(i,j) - (mv.mat_view(i,j)/mv.factor);
                }
            }

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

        // The following implementations allow to compute multiple operations.

        template <class M1, class M2>
        matrix_view<T> & operator *=(const matrix_template_addition<M1, M2> &mv)
        requires
        contain_same_element_type<typename M1::element_type, typename M2::element_type>
        and
        contain_same_element_type<T, typename M1::element_type>
        {
            if (mv.lhs.rows() != mv.rhs.rows() || mv.lhs.columns() != mv.rhs.columns())
            	throw std::domain_error("Addition operands must be the same size");

            if (ncols != mv.rhs.rows())
            	throw std::domain_error("First product operand's column number does not coincide with second operand's row number");

            matrix<T> m1{mv.lhs.rows(), mv.lhs.columns()};
            add_matrix<matrix<T>, M1, M2, T>(m1, mv.lhs, mv.rhs);
            matrix<T> m2{nrows, mv.lhs.columns()};
            mult_matrix<matrix<T>, matrix_view<T>, matrix<T>, T>(m2, *this, m1);

            *this = m2;
            return *this;
        }

        template <class M1, class M2>
        matrix_view<T> & operator *=(const matrix_template_subtraction<M1, M2> &mv)
        requires
        contain_same_element_type<typename M1::element_type, typename M2::element_type>
        and
        contain_same_element_type<T, typename M1::element_type>
        {
            if (mv.lhs.rows() != mv.rhs.rows() || mv.lhs.columns() != mv.rhs.columns())
            	throw std::domain_error("Subtraction operands must be the same size");

            if (ncols != mv.rhs.rows())
            	throw std::domain_error("First product operand's column number does not coincide with second operand's row number");

            matrix<T> m1{mv.lhs.rows(), mv.lhs.columns()};
            sub_matrix<matrix<T>, M1, M2, T>(m1, mv.lhs, mv.rhs);
            matrix<T> m2{nrows, mv.lhs.columns()};
            mult_matrix<matrix<T>, matrix_view<T>, matrix<T>, T>(m2, *this, m1);

            *this = m2;
            return *this;
        }

        template <class M>
        matrix_view<T> & operator *=(const matrix_template_scalar_product<M, T> &mv)
        requires
        contain_same_element_type<T, typename M::element_type>
        {
            if (ncols != mv.mat_view.rows())
            	throw std::domain_error("First product operand's column number does not coincide with second operand's row number");

            matrix<T> m1{mv.mat_view.rows(), mv.mat_view.columns()};
            mult_matrix_times_scalar<matrix<T>, M, T>(m1, mv.mat_view, mv.factor);
            matrix<T> m2{nrows, mv.mat_view.columns()};
            mult_matrix<matrix<T>, matrix_view<T>, matrix<T>, T>(m2, *this, m1);

            *this = m2;
            return *this;
        }

        template <class M>
        matrix_view<T> & operator *=(const matrix_template_scalar_division<M, T> &mv)
        requires
        contain_same_element_type<T, typename M::element_type>
        {
            if (ncols != mv.mat_view.rows())
            	throw std::domain_error("First product operand's column number does not coincide with second operand's row number");

            if (mv.factor == 0)
            	throw std::overflow_error("Division by 0 is not defined");

            matrix<T> m1{mv.mat_view.rows(), mv.mat_view.columns()};
            divide_matrix_by_scalar<matrix<T>, M, T>(m1, mv.mat_view, mv.factor);
            matrix<T> m2{nrows, mv.mat_view.columns()};
            mult_matrix<matrix<T>, matrix_view<T>, matrix<T>, T>(m2, *this, m1);

            *this = m2;
            return *this;
        }


        // --------------------------- SCALAR PRODUCT ---------------------------

        template <class M>
        matrix_view<T> & operator =(const matrix_template_scalar_product<M, T> &rhs)
        {
            mult_matrix_times_scalar<matrix_view<T>,M,T>(*this, rhs.mat_view, rhs.factor);
            return *this;
        }

        matrix_view<T> & operator *=(const T &rhs)
        {
            matrix_times_equal_scalar<matrix_view<T>,T>(*this, rhs);
            return *this;
        }

        // --------------------------- DIVISION ---------------------------

        template <class M>
        matrix_view<T> & operator =(const matrix_template_scalar_division<M, T> &rhs)
        {
            divide_matrix_by_scalar<matrix_view<T>,M,T>(*this, rhs.mat_view, rhs.factor);
            return *this;
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
