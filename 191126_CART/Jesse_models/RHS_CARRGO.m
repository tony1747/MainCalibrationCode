
function df = RHS_CARRGO(t, xx, p)

    
    x = xx(1); y = xx(2); 
    
    global t0_CART;
    
    tt0 = t0_CART;


    x = xx(1); y = heaviside(t-tt0)*xx(2); 
    % parameters:
    % p(1), p(2) logistic growth, capacity
    % p(3) kappa_1, killing
    % p(4) kappa_2, CART proliferation
    % p(5) theta, death rate


%     disp("inside RHS_CARRGO")
    
    
    %%%%% excluding growth parameters
    df = [ p(1)*x*(1-x/p(2)) - p(3)*x*y  ;
        p(4)*x*y - p(5)*y   ];
    


end 




