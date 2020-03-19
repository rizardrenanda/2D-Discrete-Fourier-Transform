function Iout = DCT_2D(Iin, kn)

    % input image size 
    [p, q, r] = size(Iin);
    
    % Im size divided by the kernel 
    wa = ceil(p/kn);
    ha = ceil(q/kn);
    
    % initialize output image
    Iout = zeros(p, q, r);

    % Apply and Update DCT kernel
    for k = 1 : r            
        for i = 1 : wa    
            for j = 1 : ha
                
                ws = ((i-1) * kn) + 1;
                we = min((i * kn),p);
                hs = ((j-1) * kn) + 1;
                he = min((j * kn),q);
                
                TT = Iin(ws:we,hs:he,k);
                [wk, hk] = size(TT);
                rt = zeros(wk,hk);
                
                for x1 = 1 : wk
                    for y1 = 1 : hk
                        OO = 2.0 / double(kn);
                        CR = 1;
                        CC = 1;
                        
                        % normalize
                        if x1 == 1
                            CR = 1 / sqrt(2.0);
                        end
                        % normalize
                        if y1 == 1
                            CC = 1 / sqrt(2.0);
                        end
                        
                        OO = OO * CR * CC;
                       
                        % Apply DCT for each sliding position
                        ll = 0.0;
                        for x2 = 1 : wk
                            for y2 = 1 : hk
                                temp = double(TT(x2,y2));
                                temp = temp * double(cos((((2.0*double(x2-1))+1.0)*double(x1-1)*pi)/(2.0*double(kn))));
                                temp = temp * double(cos((((2.0*double(y2-1))+1.0)*double(y1-1)*pi)/(2.0*double(kn))));
                                ll = ll + temp;
                            end
                        end
                        
                        % Store to output image
                        OO = OO * ll;
                        rt(x1,y1) = OO;
                        Iout(ws + x1-1,hs+y1-1,k) = OO;
                    end
                end
            end
        end
    end
    % Display the image in spatial and frequency domains
    figure, imshow(Iin), title('Original Image in Spatial Domain');
    figure, imshow(uint8(Iout)), title('Output Image in Frequency Domain');
end
