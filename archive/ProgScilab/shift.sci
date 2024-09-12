function us=shift(typ,u)
//{\it decale le tableau u de +1 ou -1}
//{\it avec condition aux limites periodiques}
select typ
   case '+1' then          //{\it shift +1}
           [m,n] = size(u) ;
           us(1:n-1) = u(2:n) ;
           us(n) = u(1) ;
   case '-1' then          //{\it shift -1}
           [m,n] = size(u) ;
           us(2:n) = u(1:n-1) ;
           us(1) = u(n) ;
   else 
           error('shift de +1 ou -1 ?')
end;

