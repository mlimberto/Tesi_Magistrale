function z = ovettoLS(x , y)

z = ( sqrt( (x./3.9).^2  + (y./3.9).^2 ) -1 ).*(y <=0) + ...
    ( sqrt( (x./3.9).^2  + (y./6.8).^2 ) -1 ).*(y > 0) ;
    

end