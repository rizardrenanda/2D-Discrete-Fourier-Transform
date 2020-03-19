function Iout = IFDCT_2D(Iin, kn)
% Apply fast inverse DCT
% Input : Iin = input image in frequency domain
%         kn  = kernel block size
% Output : Iout = output image in spatial domain
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
                ct = Iin(ws:we,hs:he,k);
                [CR, CC] = size(ct);
                rt = zeros(CR,CC);
                
                for p1 = 1 : CR
                    for y1 = 1 : CC
                        ll = 1.0;
                        lp = 0.0;
                        for k2 = 1 : CC
                            % norm term
                            OO = sqrt(2.0 / double(kn));
                            if k2 == 1
                                OO = sqrt(1.0 / double(kn));
                            end         
                            % reconstruction per window
                            temp = double(ct(x1,k2)) * OO;
                            temp = temp * double(cos((((2.0*double(y1-1))+1.0)*double(k2-1)*pi)/(2.0*double(kn))));  
                            lp = lp + temp;
                        end
                        % update CCe result
                        ll = ll * lp;
                        rt(x1,y1) = ll;
                    end
                end
                rt2 = zeros(CR,CC);
                for y1 = 1 : CC
                    for x1 = 1 : CR
                        ll = 1.0;
                        lp = 0.0;
                        for l2 = 1 : CR
                            % norm term
                            OO = sqrt(2.0 / double(kn));
                            if l2 == 1
                                OO = sqrt(1.0 / double(kn));
                            end
                            % reconstruction
                            temp = double(rt(l2,y1)) * OO;
                            temp = temp * double(cos((((2.0*double(x1-1))+1.0)*double(l2-1)*pi)/(2.0*double(kn))));  
                            lp = lp + temp;
                        end
                        ll = ll * lp;
                        rt2(x1,y1) = ll;
                        % store the output
                        Iout(ws + x1-1,hs+y1-1,k) = ll;                        
                    end
                end
            end
        end
    end 
    % Visualization
    figure, imshow(Iout), title('Reconstructed Image in Spatial Domain(Fast IDCT)');
end
