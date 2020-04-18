function [u_n_p_1, u_n, u_n_m_1] = Euler_step(K, M, D, h, u_n, u_n_m_1, f)
    A=((1/(h*h))*M+(1/h)*D+K);
    B=((2/(h*h))*M+(1/h)*D);
    rhs = (B*u_n-((1/(h*h))*M)*u_n_m_1+f);
    u_n_p_1 = A\rhs;
    u_n_m_1 = u_n;
    u_n     = u_n_p_1;
end