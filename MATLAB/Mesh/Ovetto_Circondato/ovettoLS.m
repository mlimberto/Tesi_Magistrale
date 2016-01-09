function z = ovettoLS(x , y)
% OVETTOLS Provides a level set function that describes the
% ovetto_circondato_sens mesh. This is used to smooth the conductivity
% coefficients

z = ( sqrt( (x./3.9).^2  + (y./3.9).^2 ) -1 ).*(y <=0) + ...
    ( sqrt( (x./3.9).^2  + (y./6.8).^2 ) -1 ).*(y > 0) ;
    

end