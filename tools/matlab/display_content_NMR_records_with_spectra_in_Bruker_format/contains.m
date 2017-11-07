function tf = contains(s, pattern)

   

        tf = findstr(s,pattern);
        if size(tf,2)>0
           tf=1;
        else
            tf=0;
        end
        return
    
end
