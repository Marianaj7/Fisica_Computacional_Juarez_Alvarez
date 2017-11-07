
__precompile__() # Este comando es para que julia precompile el paquete

module herramientas

export metodo_newton

#Documentación del método de Newton"
function metodo_newton(f,df,x0)
    z=x0
    # Se asigna x0 a z para una primera iteración
    # Se itera 100 veces
    for i in 1:100
    z=z-(f(z)/df(z));
    end
    return z
end


export rectangulo

#Documentación de los métodos de integración"
#El programa recibe tanto la función como el tamaño del intervalo, ambos antes definidos.
function rectangulo(f,b)
    p=0
    #este for hace la iteración desde 1 hasta el tamaño del intervalo -1
    for i in 1:length(b)-1
        #a corresponde al límite inferior del subintevalo, según cual sea el valor de i a elegir, de la misma manera, c corresponde al límite superior.
    a=b[i]
    c=b[i+1]
    p=p+((c-a)*f((a+c)/2))
    end
    return p
end



export trapecio
#Método del trapecio
function trapecio(f,b)
    p=0
    for i in 1:length(b)-1
    a=b[i]
    c=b[i+1]
        #esta vez se usa la fórmula del área del trapecio para la aproximación.
    p=p+((c-a)*((f(a)+f(c))/2))
    end
    return p
end



export Simpson
function Simpson(f,b)
    p=0
    for i in 1:length(b)-1
    a=b[i]
    c=b[i+1]
        #el código para calcular la integral por el método de Simpson es completamente análogo a los dos anteriores, simplemente varía en la forma de aproximarse a la integral, ahora se hace mediante un polinomio de segundo grado.
    p=p+((c-a)/6)*(f(a)+4f((a+c)/2)+f(c))
    end
    return p
end



export euler

#Documentación de los métodos de Euler"
function euler(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        x = x + f(x,t)*h
        push!(listx,x) 
     end
     return listx
end



export Euler_implicito
function metodo_Newton(f,x0,t)
    z=x0
    h=[]
    df=diff(f(x,t),x)
    for i in 1:100
    z=z-(f(z,t)/df(z,t));
        push!(h,z)
    end
    return z
    end

function Euler_implicito(f,x0,t0,t,h)
    listt=t0:h:t
    listx=[x0]
    for i in 1:length(listt)-1
         g(x)=x-x0-h*f(x,listt[i])
         x0=x0+h*f(metodo_Newton(g,Float64(x0),listt[i]),listt[i+1])
        push!(listx,x0)
    end
    return listt,listx
end
  


export runge_kutta2

#Documentación del método de Runge Kutta
function runge_kutta2(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        k1 = f(x,t);
        k2 = f(x+(h/2)*k1,t+(h/2));
        k3 = f(x+(h/2)*k2, t+(h/2));
        k4 = f(x+h*k3, t+h);
        x = x+(h/6)*(k1+2*k2+2*k3+k4);
        push!(listx,x) 
     end
     return listx
end

end
