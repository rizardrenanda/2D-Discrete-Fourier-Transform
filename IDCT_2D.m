function Iout = IDCT_2D(Iin, kn)

    % input image size 
    [p,q,r] = size(Iin);
    Iout = zeros(p,q,r);
      
    wa = ceil(p/kn);
    ha = ceil(q/kn);
    
    k = 1; i = 1; j = 1;
    for k = 1:r            
        for i= 1 : wa    
            for j = 1 : ha
                
                
                ws = ((i-1) * kn) + 1;
                we = min((i * kn),p);
                
                hs = ((j-1) * kn) + 1;
                he = min((j * kn),q);
                ct = Iin(ws:we,hs:he,k);
                [wk, hk] = size(ct);
                rt = zeros(wk,hk);
                
                
                for x1 = 1 : wk
                    for y1 = 1 : hk
                        OO = 2.0 / double(kn);

                        ll = 0.0;
                        for x2 = 1 : wk
                            for y2 = 1 : hk
                                % normalization term
                                CR = 1;
                                CC = 1;
                                if x2 == 1
                                    CR = 1.0 / sqrt(2.0);
                                end
                                if y2 == 1
                                    CC = 1.0 / sqrt(2.0);
                                end
                                
                                % Compute IDCT for each sliding position
                                temp = double(ct(x2,y2)) * CR * CC;
                                temp = temp * double(cos((((2.0*double(x1-1))+1.0)*double(x2-1)*pi)/(2.0*double(kn))));
                                temp = temp * double(cos((((2.0*double(y1-1))+1.0)*double(y2-1)*pi)/(2.0*double(kn))));
                                ll = ll + temp;
                            end
                        end
                        OO = OO * ll;
                        rt(x1,y1) = OO;
                        Iout(ws + x1-1,hs+y1-1,k) = OO;
                    end
                end
            end
        end 
    end

    figure, imshow(Iin), title('Image in Frequency Domain');
    figure, imshow(uint8(Iout)), title('Reconstructed Image in Spatial Domain');
end
