function RVNS(k_max,r_max,len_N,neighborhood_structure,e,NO_IMPROVE_LIMIT,criterio)
    r = 1;

    #Memoria vectores estaciones candidatas.
    mem_C = zeros(Int64,r_max*k_max*2,length(CANDIDATAS));
    index_mem_C = 0;

    #Solución inicial x.
    C = zeros(Int64,length(CANDIDATAS));
    E = zeros(Int64,length(ESTACIONES));
    obj = Inf;

    while true
        #println("Solución Inicial");
        C,E,obj = init_solution();
        if obj != Inf
            break;
        end
    end

    #Se guarda primera solución.
    first_C = copy(C);
    first_E = copy(E);
    first_obj = copy(obj);

    #Contador Mejoras.
    improve = 0;
    no_improve = 0;

    current_obj = zeros(Float64,r_max);
    objs_iter = [];
    cong_k = [];
    obj_k = [];
    iter_obj = [];
    objetivos_k = zeros(k_max, r_max);
    mejor_C = zeros(Int64,length(CANDIDATAS));
    mejor_E = zeros(Int64,length(ESTACIONES));

    #Comienzo iteraciones.
    while r <= r_max
        k = 1;

        println("iteración $r ====== $obj");
        push!(objs_iter,obj);

        current_obj[r] = Inf;
        mejor_obj = 999999999999999;
        mejor_obj_aux = 999999999999999;
        k_mejor = 0;

        if k_max > 1
            println("Evaluando distintos movimientos.")
            Threads.@threads for k = 1:k_max
                    #println("=== Movimiento -> k =  $(neighborhood_structure[k]) ===")

                    #Agitacion
                    #Este movimiento "shaking" devuelve una cantidad "len_N" de vecinos con c/u de los Movimientos (k)
                    array_C,array_E,array_obj = shaking(len_N,C,E,obj,neighborhood_structure[k],mem_C,index_mem_C,criterio);

                        for i=1:len_N
                            if array_obj[i] == 0 || array_obj[i] == Inf
                                break;
                            else
                                if array_obj[i] < mejor_obj_aux
                                    objetivos_k[k,r] = array_obj[i];
                                end

                                if array_obj[i] < mejor_obj
                                    mejor_C = copy(array_C[i,:]);
                                    mejor_E = copy(array_E[i,:]);
                                    mejor_obj = array_obj[i];
                                    k_mejor = k;
                                end
                            end
                        end
            end
        else
                #Agitacion
                    #Este movimiento "shaking" devuelve una cantidad "len_N" de vecinos con c/u de los Movimientos (k)
                    array_C,array_E,array_obj = shaking(len_N,C,E,obj,neighborhood_structure[k],mem_C,index_mem_C,criterio);

                    for i=1:len_N
                        if array_obj[i] == 0 || array_obj[i] == Inf
                            break;
                        else
                            if array_obj[i] < mejor_obj_aux
                                objetivos_k[k,r] = array_obj[i];
                            end

                            if array_obj[i] < mejor_obj
                                mejor_C = copy(array_C[i,:]);
                                mejor_E = copy(array_E[i,:]);
                                mejor_obj = array_obj[i];
                                k_mejor = k;
                            end
                        end
                    end
        end

        if sum(mejor_E) != 0 && mejor_obj <= current_obj[r]
            current_obj[r] = mejor_obj;
        end
        
        #Sí x' mejor que x, x = x' (Cambio de vecindario)
        if mejor_obj < obj
            println("El mejor movimiento fue con k = $(neighborhood_structure[k_mejor])")
            println("Con un valor de = $mejor_obj")
            C = copy(mejor_C);
            E = copy(mejor_E);
            obj = copy(mejor_obj);
            append!(obj_k,mejor_obj);
            push!(cong_k,neighborhood_structure[k_mejor]);
            println("=MEJORA= $obj");
            append!(iter_obj, r);
            improve += 1;
            no_improve = 0; #Si está comentado es SIN RESETEO
        else
            no_improve+=1;
            println("No hubo mejora ($no_improve)");
        end   
        
        if no_improve >= NO_IMPROVE_LIMIT
            break;
        end

        r+=1;
    end
    println("====== Resultados ======");
    println("n° iter       = $r_max");
    println("Prioridad      = $prioridad");
    println("Balance        = $balance");
    println("neighborhood_structure = $neighborhood_structure");
    println("Len_N         = $len_N");
    println("N° clusters   = $cl");
    println("N° estaciones = $(length(ESTACIONES))");
    println("N° mejoras    = $(improve)");
    println("N° iter       = $(r_max)");
    println("1° FO         = $first_obj");
    println("FO            = $obj");
    if instancia == 0
        println("Centros       = $(findall(x->x==1,C))");
    else
        ce = findall(x->x==1,C);
        centros = zeros(Int64,cl);
        x = 1;
        for i in ce
            centros[x] = ESTACIONES[i];
            x = x + 1;
        end
        println("Centros       = $centros");
    end
    println("C = $C");
    if instancia == 0
        println("E = $E");
    else
        est = zeros(Int64,length(E));
        x = 1;
        for es in E
            if es != 0
                est[x] = ESTACIONES[es];
            else
                est[x] = 0;
            end
            x = x + 1;
        end
        println("E = $est");
    end

    #Calculo de maxClusterDiameter y avgClusterDiameter
    maxClusterDiameter, avgClusterDiameter, clusterDistances = clusterDiameterStats(C, E);

    mkpath("C:/Users/joaqu/Desktop/RVNS final/Resultados finales")   
    mkpath("C:/Users/joaqu/Desktop/RVNS final/Resultados finales/$(nombre_instancia)")   
    mkpath("C:/Users/joaqu/Desktop/RVNS final/Resultados finales/$(nombre_instancia)/$(balance)_$(prioridad)")
    mkpath("C:/Users/joaqu/Desktop/RVNS final/Resultados finales/$(nombre_instancia)/$(balance)_$(prioridad)/$(r_max)/$(neighborhood_structure)_$(len_N)")     

    name = "$(balance)_$(prioridad)_exp_$(e)_$(r_max)_$(len_N)_$(improve)_$(obj)";
    filename = name*".txt"  
    open(joinpath("C:/Users/joaqu/Desktop/RVNS final/Resultados finales/$(nombre_instancia)/$(balance)_$(prioridad)/$(r_max)/$(neighborhood_structure)_$(len_N)", filename), "w") do file
        write(file, "n° iter       = $r_max \n")
        write(file, "neighborhood_structure = $neighborhood_structure \n")
        write(file, "Prioridad     = $prioridad \n")
        write(file, "balance       = $balance \n")
        write(file, "Len_N         = $len_N \n")
        write(file, "N° clusters   = $cl \n")
        write(file, "N° estaciones = $(length(ESTACIONES)) \n")
        write(file, "N° mejoras    = $(improve)\n")
        write(file, "N° iter completadas  = $(r-1)\n")
        write(file, "1° FO         = $first_obj\n")
        write(file, "FO            = $obj\n")
        if instancia == 0
            write(file, "Centros       = $(findall(x->x==1,C))\n")
        else
            ce = findall(x->x==1,C);
            println(ce)
            centros = zeros(Int64,cl);
            x = 1;
            for i in ce
                centros[x] = ESTACIONES[i];
                x = x + 1;
            end
            write(file, "Centros       = $centros\n")
        end
        write(file, "Objs= $objs_iter\n")
        write(file, "current_obj = $current_obj\n")
        write(file, "cong_K= $cong_k \n");
        write(file, "obj_k= $obj_k \n");
        write(file, "iter_obj = $iter_obj \n");
        write(file, "avgClusterDiam = $avgClusterDiameter \n");
        write(file, "maxClusterDiam = $maxClusterDiameter \n");
        write(file, "clusterDiamete = $clusterDistances   \n");
        write(file, "C = $C\n");
        if instancia == 0
            write(file, "E = $E\n");
        else
            est = zeros(Int64,length(E));
            x = 1;
            for es in E
                if es != 0
                    est[x] = ESTACIONES[es];
                else
                    est[x] = 0;
                end
                x = x + 1;
            end
            write(file, "E = $est\n");
        end
    end
    return C,E,obj,improve,obj;
end

function shaking(len_N,C,E,obj,k,mem_C,index_mem_C,criterio)
    N = zeros(Int64,len_N,length(CANDIDATAS));
    C_Arr = zeros(Int64,len_N,length(CANDIDATAS));
    E_Arr = zeros(Int64,len_N,length(CANDIDATAS));
    objArr = zeros(Int64,len_N);

    for i=1:len_N
        aux_C = zeros(Int64,length(CANDIDATAS));
        while true
            if criterio == 1
                if instancia == 0
                    aux_C = swap_center_random_grid(C,k);
                else    
                    aux_C = swap_center_random_grid(C,k);
                end
            elseif criterio == 2
                aux_C = swap_center_max_distance_grid(C,E,k);
            else
                aux_C = swap_center_priorbal_grid(C,E,k);
            end


            if compare_N(N,aux_C,len_N) && validate_connection(aux_C) && compare_N(mem_C,aux_C,index_mem_C)
                index_mem_C += 1;
                mem_C[index_mem_C,:] = aux_C;
                N[i,:] = aux_C;
                break
            end
        end
    end
    
    if len_N > 1
        Threads.@threads for i = 1:len_N
            #println("Vecino $i ====== Movimiento $k");
            #N[i,:] == aux_C#
            C_Arr[i,:] = N[i,:];
            aux_obj,aux_E = Gurobi_optimal(N[i,:]);
            if aux_obj == 0 || aux_obj == Inf
                return C_Arr,aux_E,aux_obj;
            end
            objArr[i] = aux_obj;
            E_Arr[i,:] = aux_E;
        end
    else
            i = 1;
            #println("Vecino $i ====== Movimiento $k");
            #N[i,:] == aux_C#
            C_Arr[i,:] = N[i,:];
            aux_obj,aux_E = Gurobi_optimal(N[i,:]);
            if aux_obj == 0 || aux_obj == Inf
                return C_Arr,aux_E,aux_obj;
            end
            objArr[i] = aux_obj;
            E_Arr[i,:] = aux_E;
    end
    #aux_C = N[rand(1:len_N),:];
    #aux_obj,aux_E = Gurobi_optimal(aux_C);

    #return aux_C,aux_E,aux_obj;
    return C_Arr,E_Arr,objArr;
end
