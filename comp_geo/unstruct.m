% Pseudo-unstructured mesh

% Basic mesh defined on the bi-unit square

beam_l = 16;
beam_w = 4;

% A base block has 25 nodes and 36 elements

num_blk_l = 8;
num_blk_w = 2;  

num_blk = num_blk_w*num_blk_l;

numele = 36*num_blk;

fac_w = 0.5*beam_w/num_blk_w;
fac_l = 0.5*beam_l/num_blk_l;

top = beam_w/2;

sec_1 = 18*num_blk_w + 3;

for q = 1:num_blk_l
    
    % Column 1

    for p = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(sec_1*(q-1) + 3*p-2) = 0 + 2*fac_l*(q-1);
        x(sec_1*(q-1) + 3*p-1) = 0 + 2*fac_l*(q-1);
        x(sec_1*(q-1) + 3*p)   = 0 + 2*fac_l*(q-1);
        x(sec_1*(q-1) + 3*p+1) = 0 + 2*fac_l*(q-1);
        y(sec_1*(q-1) + 3*p-2) = top - (0/3)*(fac_w) - 2*fac_w*(p-1);
        y(sec_1*(q-1) + 3*p-1) = top - (2/3)*(fac_w) - 2*fac_w*(p-1);
        y(sec_1*(q-1) + 3*p)   = top - (4/3)*(fac_w) - 2*fac_w*(p-1);
        y(sec_1*(q-1) + 3*p+1) = top - (6/3)*(fac_w) - 2*fac_w*(p-1);
        
        id(1,blocknum) = sec_1*(q-1) + 3*p-2;
        id(2,blocknum) = sec_1*(q-1) + 3*p-1;
        id(3,blocknum) = sec_1*(q-1) + 3*p;
        id(4,blocknum) = sec_1*(q-1) + 3*p+1;

    end

    % Column 2

    lenx = size(x,2);

    for p = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + 2*p-1) = (2/9)*fac_l + 2*fac_l*(q-1);
        x(lenx + 2*p)   = (2/9)*fac_l + 2*fac_l*(q-1);
        y(lenx + 2*p-1) = top - (1/3)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 2*p)   = top - (5/3)*(fac_w) - 2*fac_w*(p-1);
        
        id(5,blocknum) = lenx + 2*p-1;
        id(6,blocknum) = lenx + 2*p;
        
    end

    % Column 3

    lenx = size(x,2);

    for p = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + p) = (1/3)*fac_l + 2*fac_l*(q-1);  
        y(lenx + p) = top - fac_w - 2*fac_w*(p-1);
        
        id(7,blocknum) = lenx + p;
       
    end

    % Column 4

    lenx = size(x,2);

    for p  = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + 2*p-1) = (16/27)*fac_l + 2*fac_l*(q-1);
        x(lenx + 2*p)   = (16/27)*fac_l + 2*fac_l*(q-1);
        y(lenx + 2*p-1) = top - (16/27)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 2*p)   = top - (38/27)*(fac_w) - 2*fac_w*(p-1); 

        id(8,blocknum) = lenx + 2*p-1;
        id(9,blocknum) = lenx + 2*p;

    end

    % Column 5

    lenx = size(x,2);

    for p  = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + p)     = (2/3)*fac_l + 2*fac_l*(q-1);
        x(lenx + p+1)   = (2/3)*fac_l + 2*fac_l*(q-1);
        y(lenx + p)     = top - (0)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + p+1)   = top - (2)*(fac_w) - 2*fac_w*(p-1);

        id(10,blocknum) = lenx + p;
        id(11,blocknum) = lenx + p+1;

    end

    % Column 6

    lenx = size(x,2);

    for p  = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + 3*p-2) = (0.9)*fac_l + 2*fac_l*(q-1);
        x(lenx + 3*p-1) = (1)*fac_l + 2*fac_l*(q-1);
        x(lenx + 3*p)   = (1.1)*fac_l + 2*fac_l*(q-1);
        y(lenx + 3*p-2) = top - (1/3)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 3*p-1) = top - (1)*(fac_w)   - 2*fac_w*(p-1);
        y(lenx + 3*p)   = top - (5/3)*(fac_w) - 2*fac_w*(p-1); 

        id(12,blocknum) = lenx + 3*p-2;
        id(13,blocknum) = lenx + 3*p-1;
        id(14,blocknum) = lenx + 3*p;

    end
   
    % Column 7

    lenx = size(x,2);

    for p  = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + p)     = (4/3)*fac_l + 2*fac_l*(q-1);
        x(lenx + p+1)   = (4/3)*fac_l + 2*fac_l*(q-1);
        y(lenx + p)     = top - (0)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + p+1)   = top - (2)*(fac_w) - 2*fac_w*(p-1); 

        id(15,blocknum) = lenx + p;
        id(16,blocknum) = lenx + p+1;

    end
    
    % Column 8

    lenx = size(x,2);

    for p  = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + 2*p-1) = (38/27)*fac_l + 2*fac_l*(q-1);
        x(lenx + 2*p)   = (38/27)*fac_l + 2*fac_l*(q-1);
        y(lenx + 2*p-1) = top - (16/27)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 2*p)   = top - (38/27)*(fac_w) - 2*fac_w*(p-1); 

        id(17,blocknum) = lenx + 2*p-1;
        id(18,blocknum) = lenx + 2*p;

    end

    % Column 9

    lenx = size(x,2);

    for p = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + p) = (5/3)*fac_l + 2*fac_l*(q-1);  
        y(lenx + p) = top - fac_w - 2*fac_w*(p-1);

        id(19,blocknum) = lenx + p;
       
    end

    % Column 10

    lenx = size(x,2);

    for p = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + 2*p-1) = (16/9)*fac_l + 2*fac_l*(q-1);
        x(lenx + 2*p)   = (16/9)*fac_l + 2*fac_l*(q-1);
        y(lenx + 2*p-1) = top - (1/3)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 2*p)   = top - (5/3)*(fac_w) - 2*fac_w*(p-1); 

        id(20,blocknum) = lenx + 2*p-1;
        id(21,blocknum) = lenx + 2*p;
        
    end
    
    % Column 11
    
    lenx = size(x,2);

    for p = 1:num_blk_w
        
        blocknum = num_blk_w*(q-1) + p;
    
        x(lenx + 3*p-2) = 2*fac_l + 2*fac_l*(q-1);
        x(lenx + 3*p-1) = 2*fac_l + 2*fac_l*(q-1);
        x(lenx + 3*p)   = 2*fac_l + 2*fac_l*(q-1);
        x(lenx + 3*p+1) = 2*fac_l + 2*fac_l*(q-1);
        y(lenx + 3*p-2) = top - (0/3)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 3*p-1) = top - (2/3)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 3*p)   = top - (4/3)*(fac_w) - 2*fac_w*(p-1);
        y(lenx + 3*p+1) = top - (6/3)*(fac_w) - 2*fac_w*(p-1);
        
        id(22,blocknum) = lenx + 3*p-2;
        id(23,blocknum) = lenx + 3*p-1;
        id(24,blocknum) = lenx + 3*p;
        id(25,blocknum) = lenx + 3*p+1;
    
    end

end

% Base block connectivity array

nodet = [];

for i = 1:num_blk

nodet = [ nodet;
          id(1,i)  id(5,i)  id(10,i);
          id(2,i)  id(5,i)  id(1,i);
          id(5,i)  id(2,i)  id(8,i);
          id(5,i)  id(8,i)  id(10,i);
          id(10,i) id(8,i)  id(12,i);
          id(10,i) id(12,i) id(15,i);
          id(12,i) id(17,i) id(15,i);
          id(15,i) id(17,i) id(20,i);
          id(15,i) id(20,i) id(22,i);
          id(20,i) id(23,i) id(22,i);
          id(20,i) id(17,i) id(23,i);
          id(17,i) id(19,i) id(23,i);
          id(19,i) id(24,i) id(23,i);
          id(19,i) id(18,i) id(24,i);
          id(18,i) id(21,i) id(24,i);
          id(21,i) id(25,i) id(24,i);
          id(16,i) id(25,i) id(21,i);
          id(18,i) id(16,i) id(21,i);
          id(14,i) id(16,i) id(18,i);
          id(11,i) id(16,i) id(14,i);
          id(9,i)  id(11,i) id(14,i);
          id(6,i)  id(11,i) id(9,i);
          id(4,i)  id(11,i) id(6,i);
          id(3,i)  id(4,i)  id(6,i);
          id(3,i)  id(6,i)  id(9,i);
          id(3,i)  id(9,i)  id(7,i);
          id(3,i)  id(7,i)  id(2,i);
          id(2,i)  id(7,i)  id(8,i);
          id(7,i)  id(13,i) id(8,i);
          id(7,i)  id(9,i)  id(13,i);
          id(9,i)  id(14,i) id(13,i);
          id(13,i) id(14,i) id(18,i);
          id(19,i) id(13,i) id(18,i);
          id(19,i) id(17,i) id(13,i);
          id(13,i) id(17,i) id(12,i);
          id(13,i) id(12,i) id(8,i)];
                   
end
          
node = nodet';
          
numnod = size(x,2);

clear nodet num_blk_l num_blk_w num_blk blocknum p q i lenx id
clear fac_l fac_w sec_1 top beam_l beam_w