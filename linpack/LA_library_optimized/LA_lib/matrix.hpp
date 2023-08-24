//==================================================================================================
//  File:       matrix.hpp
//
//  Summary:    This header defines the matrix class template.
//==================================================================================================

#ifndef LINEAR_ALGEBRA_MATRIX_HPP
#define LINEAR_ALGEBRA_MATRIX_HPP

namespace LA_lib {

    template <class M1, class M2>
    struct matrix_template_addition {
        const M1 & lhs;
        const M2 & rhs;
    };

    template <class M1, class M2>
    struct matrix_template_subtraction {
        const M1 & lhs;
        const M2 & rhs;
    };

    template <class M, typename T>
    struct matrix_template_scalar_product {
        const M & mat_view;
        T factor;
    };

    template <class M, typename T>
    struct matrix_template_scalar_division {
        const M & mat_view;
        T factor;
    };


    template <typename T1, typename T2>
    concept contain_same_element_type = std::is_same<T1, T2>::value;

    template <class M1, class M2>
    matrix_template_addition<M1, M2> operator +(const M1 & lhs, const M2 & rhs)
    requires
    contain_same_element_type<typename M1::element_type, typename M2::element_type>
    {
        return matrix_template_addition<M1, M2>{lhs,rhs};
    }

    template <class M1, class M2>
    matrix_template_subtraction<M1, M2> operator -(const M1 & lhs, const M2 & rhs)
    requires
    contain_same_element_type<typename M1::element_type, typename M2::element_type>
    {
        return matrix_template_subtraction<M1, M2>{lhs,rhs};
    }

    template <class M, typename T>
    matrix_template_scalar_product<M,T> operator *(const M & m, T n)
    requires
    contain_same_element_type<typename M::element_type, T>
    {
        return matrix_template_scalar_product<M,T>{m,n};
    }

    template <class M, typename T>
    matrix_template_scalar_product<M,T> operator *(T n, const M & m)
    requires
    contain_same_element_type<typename M::element_type, T>
    {
        return matrix_template_scalar_product<M,T>{m,n};
    }

    template <class M, typename T>
    matrix_template_scalar_division<M,T> operator /(const M & m, T n)
    requires
    contain_same_element_type<typename M::element_type, T>
    {
        return matrix_template_scalar_division<M, T>{m,n};
    }

    template <typename T>
    class matrix_view;

    template <class M>
    class matrix_transpose;

    template<typename T>
    class matrix {
    private:

        std::vector <T> buffer;
        size_t nrows;
        size_t ncols;

    public:

        using element_type = T;

// ==============================================================================
//                         MATRIX OBJECT INITIALIZATION
// ==============================================================================

        matrix() :
                buffer(0),
                nrows{0},
                ncols{0} {}

        template<typename N>
        matrix(N rows, N cols):
                buffer(rows * cols),
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}

        template<typename N>
        matrix(std::initializer_list <T> &rhs, N rows, N cols):
                buffer{rhs},
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}

        template<typename N>
        matrix(std::initializer_list <T> &&rhs, N rows, N cols):
                buffer{rhs},
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}

        template<typename N>
        matrix(std::vector <T> &rhs, N rows, N cols):
                buffer{rhs},
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}


// ==============================================================================
//                              MATRIX DATA ACCESS
// ==============================================================================

        T operator()(int ri, int ci) const
        {
            return buffer[(ri * ncols) + ci];
        }

        T &operator()(int ri, int ci)
        {
            return buffer[(ri * ncols) + ci];
        }

        T operator()(int i) const
        {
            if (nrows > 1 && ncols > 1)
                throw std::invalid_argument("2D object cannot be accessed by single index");

            return buffer[i];
        }

        T &operator()(int i)
        {
            if (nrows > 1 && ncols > 1)
                throw std::invalid_argument("2D object cannot be accessed by single index");

            return buffer[i];
        }

        T *data() { return buffer.data(); }


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
        size() const noexcept { return (buffer.size()); }

        constexpr void
        resize(size_t rn, size_t cn)
        {
            if (rn == 0 && cn != 0)
                throw std::invalid_argument("Matrix cannot have zero rows and non-zero column number");
            if (rn != 0 && cn == 0)
                throw std::invalid_argument("Matrix cannot have zero rows and non-zero column number");

            if (rn != nrows || cn != ncols) {
                buffer.reserve(rn * cn);
                buffer.resize(rn * cn);
                buffer.shrink_to_fit();
                nrows = rn;
                ncols = cn;
            }
        }


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


// ==============================================================================
//                         ASSIGNMENT OPERATOR OVERLOADING
// ==============================================================================

        matrix<T> &operator=(matrix<T> const &rhs)
        {
            if (this != &rhs) {
                this->resize(rhs.nrows, rhs.ncols);
                buffer = rhs.buffer;
                nrows = rhs.nrows;
                ncols = rhs.ncols;
                return *this;
            } else {
                return *this;
            }
        }

        matrix<T> &operator=(std::initializer_list <T> const &rhs)
        {
            if (!std::empty(rhs)) {
                this->resize(1, rhs.size());
                buffer = rhs;
                nrows = 1;
                ncols = rhs.size();
                return *this;
            } else {
                this->resize(0, 0);
                buffer = rhs;
                nrows = 0;
                ncols = 0;
                return *this;
            }

        }

        matrix<T> &operator=(std::initializer_list <std::initializer_list<T>> const &rhs)
        {
            if (!check_equal_size_rows(rhs))
                throw std::domain_error("The number of columns must be equal in every row");

            if (is_list_empty(rhs)) {
                this->resize(0, 0);
                buffer = {};
                nrows = 0;
                ncols = 0;
                return *this;
            } else {
                this->resize(rhs.size(), rhs.begin()->size());

                int i = 0;
                for (auto const &v1: rhs) {
                    for (auto const &v2: v1) {
                        buffer[i] = v2;
                        ++i;
                    }
                }

                nrows = rhs.size();
                ncols = rhs.begin()->size();
                return *this;
            }
        }

        matrix<T> &operator=(matrix_view<T> const &rhs)
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
        matrix<T> &operator=(matrix_transpose<M> const &rhs)
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
        matrix<T> & operator =(const matrix_template_addition<M1, M2> &mv)
        requires
                contain_same_element_type<T, typename M1::element_type>
        {
            add_matrix<matrix<T>, M1, M2, T>(*this, mv.lhs, mv.rhs);
            return *this;
        }

        template <class M>
        matrix<T> & operator +=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            plus_equal_matrix<matrix<T>, M>(*this, rhs);
            return *this;
        }

        // The following implementations allow to compute multiple operations in an efficient way.
        // No temporal objects are generated to obtain the result.

        template <class M1, class M2>
        matrix<T> & operator +=(const matrix_template_addition<M1, M2> &mv)
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
        matrix<T> & operator +=(const matrix_template_subtraction<M1, M2> &mv)
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
        matrix<T> & operator +=(const matrix_template_scalar_product<M, T> &mv)
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
        matrix<T> & operator +=(const matrix_template_scalar_division<M, T> &mv)
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
        matrix<T> & operator =(const matrix_template_subtraction<M1, M2> &mv)
        requires
                contain_same_element_type<T, typename M1::element_type>
        {
            sub_matrix<matrix<T>,M1, M2, T>(*this, mv.lhs, mv.rhs);
            return *this;
        }

        template <class M>
        matrix<T> & operator -=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            minus_equal_matrix<matrix<T>, M>(*this, rhs);
            return *this;
        }

        // The following implementations allow to compute multiple operations in an efficient way.
        // No temporal objects are generated to obtain the result.

        template <class M1, class M2>
        matrix<T> & operator -=(const matrix_template_addition<M1, M2> &mv)
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
        matrix<T> & operator -=(const matrix_template_subtraction<M1, M2> &mv)
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
        matrix<T> & operator -=(const matrix_template_scalar_product<M, T> &mv)
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
        matrix<T> & operator -=(const matrix_template_scalar_division<M, T> &mv)
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
        {
            matrix<T> m;
            mult_matrix<matrix<T>, M1, M2, T>(m, lhs, rhs);
            return m;
        }

        template <class M>
        matrix<T> operator *=(const M & rhs)
        requires
        contain_same_element_type<typename M::element_type, T>
        {
            matrix<T> m;
            mult_matrix<matrix<T>, matrix<T>, M, T>(m, *this, rhs);
            return m;
        }

        // The following implementations allow to compute multiple operations.

        template <class M1, class M2>
        matrix<T> & operator *=(const matrix_template_addition<M1, M2> &mv)
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

            return m2;
        }

        template <class M1, class M2>
        matrix<T> & operator *=(const matrix_template_subtraction<M1, M2> &mv)
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

            return m2;
        }

        template <class M>
        matrix<T> & operator *=(const matrix_template_scalar_product<M, T> &mv)
        requires
        contain_same_element_type<T, typename M::element_type>
        {
            if (ncols != mv.mat_view.rows())
                throw std::domain_error("First product operand's column number does not coincide with second operand's row number");

            matrix<T> m1{mv.mat_view.rows(), mv.mat_view.columns()};
            mult_matrix_times_scalar<matrix<T>, M, T>(m1, mv.mat_view, mv.factor);
            matrix<T> m2{nrows, mv.mat_view.columns()};
            mult_matrix<matrix<T>, matrix_view<T>, matrix<T>, T>(m2, *this, m1);

            return m2;
        }

        template <class M>
        matrix<T> & operator *=(const matrix_template_scalar_division<M, T> &mv)
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

            return m2;
        }
        

        // --------------------------- SCALAR PRODUCT ---------------------------

        template <class M>
        matrix<T> & operator =(const matrix_template_scalar_product<M, T> &rhs)
        {
            mult_matrix_times_scalar<matrix<T>,M,T>(*this, rhs.mat_view, rhs.factor);
            return *this;
        }

        matrix<T> & operator *=(const T &rhs)
        {
            matrix_times_equal_scalar<matrix<T>,T>(*this, rhs);
            return *this;
        }

        // --------------------------- DIVISION ---------------------------

        template <class M>
        matrix<T> & operator =(const matrix_template_scalar_division<M, T> &rhs)
        {
            divide_matrix_by_scalar<matrix<T>,M,T>(*this, rhs.mat_view, rhs.factor);
            return *this;
        }

        matrix<T> & operator /=(const T &rhs)
        {
            matrix_divide_equal_scalar<matrix<T>,T>(*this, rhs);
            return *this;
        }



        // --------------------------- COMPARISON ---------------------------

        matrix<T> operator==(matrix<T> const &rhs)
        {
            return (this->buffer == rhs.buffer);
        }

        matrix<T> operator!=(matrix<T> const &rhs)
        {
            return (this->buffer != rhs.buffer);
        }


// ==============================================================================
//                               SUBELEMENT ACCESS
// ==============================================================================

        matrix_view<T> submatrix(int ri, int rn, int ci, int cn)
        {
            return matrix_view(*this, ri, rn, ci, cn);
        }

        matrix_view<T> column(int i)
        {
            return matrix_view(*this, 0, nrows, i, 1);
        }

        matrix_view<T> row(int i)
        {
            return matrix_view(*this, i, 1, 0, ncols);
        }

        matrix_transpose<matrix<T>> t(){
            return matrix_transpose<matrix<T>>{*this};
        }
    };

}

#endif //LINEAR_ALGEBRA_MATRIX_HPP
