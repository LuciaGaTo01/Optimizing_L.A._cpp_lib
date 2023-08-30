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

    template <typename T>
    class matrix_view;

    template <class M>
    class matrix_transpose;

    template <typename T1, typename T2>
    concept contain_same_element_type = std::is_same<T1, T2>::value;

    template <class M1, class M2>
    class matrix_template_addition {
    private:
        const M1 &lhs;
        const M2 &rhs;
    public:
        using element_type = typename M1::element_type;

        matrix_template_addition(const M1 &lhs, const M2 &rhs):
                lhs{lhs},
                rhs{rhs}
                {assert((lhs.rows() == rhs.rows() && lhs.columns() == rhs.columns()
                         && "Both matrices must have the same size"));}

        [[nodiscard]] constexpr size_t
        rows() const noexcept { return lhs.rows(); }

        [[nodiscard]] constexpr size_t
        columns() const noexcept { return lhs.columns(); }

        [[nodiscard]] constexpr size_t
        size() const noexcept { return lhs.size(); }

        element_type operator()(int ri, int ci) const
        {
            return lhs(ri,ci) + rhs(ri,ci);
        }
    };

    template <class M1, class M2>
    class matrix_template_subtraction {
    private:
        const M1 &lhs;
        const M2 &rhs;
    public:
        using element_type = typename M1::element_type;

        matrix_template_subtraction(const M1 &lhs, const M2 &rhs):
                lhs{lhs},
                rhs{rhs}
                {assert((lhs.rows() == rhs.rows() && lhs.columns() == rhs.columns()
                         && "Both matrices must have the same size"));}

        [[nodiscard]] constexpr size_t
        rows() const noexcept { return lhs.rows(); }

        [[nodiscard]] constexpr size_t
        columns() const noexcept { return lhs.columns(); }

        [[nodiscard]] constexpr size_t
        size() const noexcept { return lhs.size(); }

        element_type operator()(int ri, int ci) const
        {
            return lhs(ri,ci) - rhs(ri,ci);
        }
    };

    template <class M1, class M2>
    class matrix_template_product {
    private:
        const M1 &lhs;
        const M2 &rhs;
    public:
        using element_type = typename M1::element_type;

        matrix_template_product(const M1 &lhs, const M2 &rhs):
                lhs{lhs},
                rhs{rhs}
                {assert((lhs.columns() == rhs.rows()
                         && "First matrix columns must be equal to second matrix rows"));}

        [[nodiscard]] constexpr size_t
        rows() const noexcept { return lhs.rows(); }

        [[nodiscard]] constexpr size_t
        columns() const noexcept { return rhs.columns(); }

        [[nodiscard]] constexpr size_t
        size() const noexcept { return lhs.rows() * rhs.columns(); }

        element_type operator()(int ri, int ci) const
        {
            element_type temp = 0;
            for (int i = 0; i < lhs.columns(); ++i)
                temp += lhs(ri,i) * rhs(i,ci);
            return temp;
        }
    };

    template <class M1, typename T>
    class matrix_template_scalar_product {
    private:
        const M1 &lhs;
        const T rhs;
    public:
        using element_type = T;

        matrix_template_scalar_product(const M1 &lhs, const T rhs):
                lhs{lhs},
                rhs{rhs}{}

        [[nodiscard]] constexpr size_t
        rows() const noexcept { return lhs.rows(); }

        [[nodiscard]] constexpr size_t
        columns() const noexcept { return lhs.columns(); }

        [[nodiscard]] constexpr size_t
        size() const noexcept { return lhs.size(); }

        element_type operator()(int ri, int ci) const
        {
            return lhs(ri,ci) * rhs;
        }
    };

    template <class M1, typename T>
    class matrix_template_scalar_division {
    private:
        const M1 &lhs;
        const T rhs;
    public:
        using element_type = T;

        matrix_template_scalar_division(const M1 &lhs, const T rhs):
                lhs{lhs},
                rhs{rhs}{}

        [[nodiscard]] constexpr size_t
        rows() const noexcept { return lhs.rows(); }

        [[nodiscard]] constexpr size_t
        columns() const noexcept { return lhs.columns(); }

        [[nodiscard]] constexpr size_t
        size() const noexcept { return lhs.size(); }

        element_type operator()(int ri, int ci) const
        {
            return lhs(ri,ci)/rhs;
        }
    };

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

    template <class M1, class M2>
    matrix_template_product<M1,M2> operator *(const M1 & lhs, const M2 &rhs)
    requires
    contain_same_element_type<typename M1::element_type, typename M2::element_type>
    {
        return matrix_template_product<M1,M2>{lhs,rhs};
    }

    template <class M, typename T>
    matrix_template_scalar_product<M,T> operator *(const M & m, const T n)
    requires
    contain_same_element_type<typename M::element_type, T>
    {
        return matrix_template_scalar_product<M,T>{m,n};
    }

    template <class M, typename T>
    matrix_template_scalar_product<M,T> operator *(const T n, const M & m)
    requires
    contain_same_element_type<typename M::element_type, T>
    {
        return matrix_template_scalar_product<M,T>{m,n};
    }

    template <class M, typename T>
    matrix_template_scalar_division<M,T> operator /(const M & m, const T n)
    requires
    contain_same_element_type<typename M::element_type, T>
    {
        return matrix_template_scalar_division<M, T>{m,n};
    }

    template <class M1, class M2, class M3, typename T>
    void mult_matrix (M1 &m, M2 const &lhs, M3 const &rhs)
    {
        assert((lhs.columns() == rhs.rows()
                && "First matrix column number must be equal to second matrix row number"));

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
        assert((nrows == rhs.rows() && ncols == rhs.columns()
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
