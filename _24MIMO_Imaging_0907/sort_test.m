function [x1,y1,inds] = sort_test(x,y,x0,y0)
    xx = x;
    yy = y;
    I = size(x0,2);
    x1 = zeros(1,I);
%     x1_tmp = -inf*ones(K,I);
    y1 = zeros(1,I);
    inds_ = [];
%     y1_tmp = -inf*ones(K,I);
    for i = 1:I
        R = rangeangle([x;y;zeros(size(x))],[x0(i);y0(i);0]);
        if min(R) < 1
            ind = find(R==min(R),1);
            x1(i) = x(ind);
            y1(i) = y(ind); %#ok<*FNDSB>
            x(ind) = [];
            y(ind) = [];
%             if ind == 1
%                 x = x(2:end);
%                 y = y(2:end);
%             else if ind == I %#ok<SEPEX>
%                     x = x(1:end-1);
%                     y = y(1:end-1);
%                 else 
%                     x = [x(1:(ind-1)) x((ind+1):end)]; %#ok<NASGU>
%                     x = [y(1:(ind-1)) y((ind+1):end)];
%                 end
%             end
        else
            inds_ = [inds_,i];
        end
    end
    for i = 1:length(x)
         R = rangeangle([x0(inds_);y0(inds_);zeros(1,length(inds_))],[x(i);y(i);0]);
         ind2 = find(R==min(R),1);
         x1(inds_(ind2)) = x(i);
         y1(inds_(ind2)) = y(i);
         inds_(ind2) = [];
    end
    for i = 1:I
        inds(i) = find(x1==xx(i),1);
    end
end
                
                
    
            