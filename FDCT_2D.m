function Iout = FDCT_2D(Iin, kn)
% Apply fast DCT
% Input : Iin = input image in SPATIAL domain
%         kn  = kernel block size
% Output : Iout = output image in frequency domain
% input image size 
[p, q, r] = size(Iin);
Iout = zeros(p, q, r);
% Im size divided by the kernel 
wa = ceil(p/kn);
ha = ceil(q/kn);
    
    for k = 1:r            
        for i = 1 : wa   
            for j = 1 : ha
                ws = ((i-1) * kn) + 1;
                we = min((i * kn),p);
                
                hs = ((j-1) * kn) + 1;
                he = min((j * kn),q);
                % padding
                ct = I(ws:we,hs:he,k);
                [CC, CR] = size(ct);
                rt = zeros(CC,CR);

                for x1 = 1 : CC
                    for y1 = 1 : CR
                        % norm term
                        OO = sqrt(2.0 / double(kn));
                        if y1 == 1
                            OO = sqrt(1.0 / double(kn));
                        end
                        ll = OO;
                        lp = 0.0;
                        % compute DCT
                        for k2 = 1 : CR
                            temp = double(ct(x1,k2));
                            temp = temp * double(cos((((2.0*double(k2-1))+1.0)*double(y1-1)*pi)/(2.0*double(kn))));  
                            lp = lp + temp;
                        end
                        ll = ll * lp;
                        rt(x1,y1) = ll;
                    end
                end

                rt2 = zeros(CC,CR);
                for y1 = 1 : CR
                    for x1 = 1 : CC
                        % norm term
                        OO = sqrt(2.0 / double(kn));
                        if y1 == 1
                            OO = sqrt(1.0 / double(kn));
                        end
                        ll = OO;
                        lp = 0.0;
                        % compute DCT
                        for l2 = 1 : CC
                            temp = double(rt(l2,y1));
                            temp = temp * double(cos((((2.0*double(l2-1))+1.0)*double(x1-1)*pi)/(2.0*double(kn))));  
                            lp = lp + temp;
                        end
                        ll = ll * lp;
                        rt2(x1,y1) = ll;
                        % store output
                        Iout(ws + x1-1,hs+y1-1,k) = ll;                        
                    end
                end
            end
        end
    end
    % Visualization
    figure, imshow(Iout), title('Output Image in Frequency Domain (Fast DCT)');
end
