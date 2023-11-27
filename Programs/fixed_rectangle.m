% Produces random rectangles on a square array of specified size.
function [arrayvector,arrayvector2] = fixed_rectangle(n)
    array = ones(n+1);
    array2 = zeros(n+1);
    x = linspace(0,1,n+1);
    y = x;
    [xx,yy] = meshgrid(x,y);
    number_rect = 15;
    randomnums = importdata('Random_Numbers_Initial_Rectangles.mat');
    randomnumcount = number_rect*3;
    for i=1:number_rect
        upperx = randomnums(3*(i-1)+1);
        lowerx = upperx^2;
        uppery = randomnums(3*(i-1)+2);
        lowery = uppery^2;
        value = randomnums(3*(i-1)+3);
        for i1=1:n+1
            for i2=1:n+1
                if (xx(i1,i2)>=lowerx) && (xx(i1,i2)<=upperx) && (yy(i1,i2)>=lowery) && (yy(i1,i2)<=uppery)
                    array(i1,i2) = value;
                end
            end
        end
    end
    array = array(2:end,2:end);
    arrayvector = reshape(array,n^2,1);
    %arrayvector2 = ones(n^2,1);
    
    for i=1:number_rect
        upperx = randomnums(randomnumcount + 3*(i-1)+1);
        lowerx = upperx^2;
        uppery = randomnums(randomnumcount + 3*(i-1)+2);
        lowery = uppery^2;
        value = randomnums(randomnumcount + 3*(i-1)+3);
        for i1=1:n+1
            for i2=1:n+1
                if (xx(i1,i2)>=lowerx) && (xx(i1,i2)<=upperx) && (yy(i1,i2)>=lowery) && (yy(i1,i2)<=uppery)
                    array2(i1,i2) = value;
                end
            end
        end
    end
    array2 = array2(2:end,2:end);
    arrayvector2 = reshape(array2,n^2,1);
    %save('Random_Rectangles_n50.mat','arrayvector','arrayvector2')
end