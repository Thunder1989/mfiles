function foo(var)
    a = struct();
    a.x = var;
    a = bar(a);
    a.x
    a = bar(a);
    a.x

function a = bar(a)
    a.x = a.x + 1 ;
   