function [x,y] = tadpolegeo(M)
    H_tail = 1/8;
    R_head = 1/8;
    L = 3/4;
    C_head_x =  1 - (1 -L)/2 - R_head;
    C_head_y = 1/2;

    X_off = (R_head^2-(H_tail/2)^2)^(1/2);
    theta_off = atan(H_tail/(2*X_off))
	L_side = C_head_x - X_off -(1 - L)/2;
    L_head = 2*(pi - theta_off)*R_head;
    perim = H_tail + L_head + 2*L_side;


    M_bot = floor(M*(L_side/perim));
	M_head = floor(M*(L_head/perim));
	M_top = M_bot;
	M_tail = M-M_bot-M_top-M_head;
    
    ds_y = H_tail/M_tail;
	ds_x = 0;
	x(1) = (1. - L)/2.
	y(1) = (1. - H_tail + ds_y)/2.
	
    
    for i = 2:M_tail
        x(i) = x(i-1) + ds_x;
		y(i) = y(i-1) + ds_y;
        
    end
    
    
    ds_x = L_side/M_top;
	ds_y = 0;
	x(M_tail+1) = (1. - L + ds_x)/2;
	y(M_tail+1) = (1. + H_tail + ds_y)/2;
    
    for I = M_tail+2:M_tail+M_top
		x(I) = x(I-1) + ds_x;
		y(I) = y(I-1) + ds_y;
    end
    
    d_theta = 2*(pi - theta_off)/M_head;
  	
    for I = 1:M_head
		K = I+M_tail+M_top
		theta = pi - theta_off - (I-1./2)*d_theta;
		x(K) = C_head_x + R_head*cos(theta);
		y(K) = C_head_y + R_head*sin(theta);
    end
    	
    ds_x = -L_side/M_bot;
	ds_y = 0;
	x(M_tail+M_top+M_head+1) = C_head_x - X_off + ds_x/2.
	y(M_tail+M_top+M_head+1) = (1. - H_tail + ds_y)/2.
    
    for I = M_tail+M_top+M_head+2:M
		x(I) = x(I-1) + ds_x;
		y(I) = y(I-1) + ds_y;
    end
end