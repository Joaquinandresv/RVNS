
function init_solution()
    _C = zeros(Int64,length(CANDIDATAS));
    _E = zeros(Int64,length(ESTACIONES));
    _obj = Inf;

    #TESTING
    _Et = zeros(Int64,length(ESTACIONES));
    _objt = Inf;

    _Ett = zeros(Int64,length(ESTACIONES));
    _objtt = Inf;

    if instancia == 0
        while true
            _C = init_solution_C_grid();
            if validate_connection(_C)
                break;
            end
        end
    else
        while true
            _C = init_solution_C_grid();
            if validate_connection(_C)
                break;
            end
        end
    end
    _obj,_E = Gurobi_optimal(_C); #@time Gurobi_optimal(_C);

    return _C,_E,_obj;
end

function init_solution_C()
    println("Solución inicial sin grilla");
    _C = zeros(Int64,length(CANDIDATAS));
    for i=1:cl
        while true
            y = rand(CANDIDATAS);
            index = findall(x->x==y, CANDIDATAS)[1]
            if (~(_C[index] == 1) && ~(length(findall(x->x==1,_C)) > cl))
                _C[index] = 1;
                break;
            end
        end
    end
    return _C;
end

function init_solution_C_grid()
    println("Solución inicial con grilla");
    _C = zeros(Int64,length(CANDIDATAS));
    for i=1:cl
        zone = findall(x->trunc(x)==1,vec(M[i,:])); ##RETORNA LA ZONA, Parte con i=1; saca la zona
            if !isempty(zone)                               ## y sus respectivas estaciones, "==1"
                while true  
                    x = rand(zone);
                    if (~(_C[x] == 1) && ~(length(findall(x->x==1,_C)) > cl))
                        _C[x] = 1;
                        break;
                    end
                end
            end
    end
    return _C;
end

function compare_N(N,C,len_N)
    cond = 0;
    for i=1:len_N
        if N[i,:] == C[:]
            cond = 1;
        end
    end
    if (cond == 0)
        return true;
    else
        return false;
    end
end

function init_prior()
    _P = zeros(Int64,length(ESTACIONES));
    for i in ESTACIONES
        prior_aux = PRIORIDADES[findall(x->x==1,vec(prior[i,:]'))];
        _P[i] = prior_aux;
    end
    return _P;
end

function swap_center_random(C,k)
    _C = zeros(Int64,length(CANDIDATAS));
    _C = copy(C);
    centers_old = findall(x->x==1,C)
    a = [];
    b = [];
    for i=1:k
        centers = (x->x==1,_C);
        e = 0;
        while true
            e = rand(centers_old);
            if (~(e in a))
                append!(a,e);
                break;
            end
        end
        index=0;
        while true
            y = rand(CANDIDATAS);
            index = findall(x->x==y, CANDIDATAS)[1]
            if (~(index in centers) && ~(index in b) && ~(index in centers_old))
                append!(b,index);
                break;
            end
        end
        _C[a[i]] = 0;
        _C[index] = 1;
    end
    return _C;
end

function swap_center_random_grid(C,k)
    _C = zeros(Int64,length(CANDIDATAS));
    _C = copy(C);
    centers_old = findall(x->x==1,C)
    a = [];
    b = [];

    for i=1:k
            zone = 0;
            centers = (x->x==1,_C);
            e = 0;
            while true
                e = rand(centers_old);
                if (~(e in a))
                    append!(a,e);
                    zone = get_zone(e);
                    break;
                end
            end
            x=0;
            while true
                x = rand(findall(x->trunc(x)==1,vec(M[zone,:])));
                if (~(x in centers) && ~(x in b) && ~(x in centers_old))
                    append!(b,x);
                    break;
                end
            end
            _C[a[i]] = 0;
            _C[x] = 1;
    end
    return _C;
end


function swap_center_max_distance(C,E,k)
    _C = zeros(Int64,length(CANDIDATAS));
    _C = copy(C);
    centers_old = findall(x->x==1,C)
    a = [];
    b = [];
    for i=1:k
        centers = (x->x==1,_C);
        e = 0;
        matrix_distances = get_distances_cluster_system(C,E)
        e = trunc(Int,matrix_distances[i,1]);
        append!(a,e);

        x=0;
        x = rand(findall(x-> ~(x in b) && ~(x in centers_old) && ~(x in centers),CANDIDATAS))
        append!(b,x);

        _C[e] = 0;
        _C[x] = 1;
    end
    return _C;
end

function swap_center_max_distance_grid(C,E,k)
    _C = zeros(Int64,length(CANDIDATAS));
    _C = copy(C);
    centers_old = findall(x->x==1,C)
    a = [];
    b = [];
    for i=1:k
        centers = (x->x==1,_C);
        e = 0;
        matrix_distances = get_distances_cluster_system(C,E)
        e = trunc(Int,matrix_distances[i,1]);
        append!(a,e);
        zone =get_zone(e);

        x=0;
        while true
            x = rand(findall(x->trunc(x)==1,vec(M[zone,:])));
            if (~(x in centers) && ~(x in b) && ~(x in centers_old))
                append!(b,x);
                break;
            end
        end

        _C[e] = 0;
        _C[x] = 1;
    end
    return _C;
end


function swap_center_priorbal(C,E,k)
    _C = zeros(Int64,length(CANDIDATAS));
    _C = copy(C);
    centers_old = findall(x->x==1,C)
    a = [];
    b = [];
    for i=1:k
        centers = (x->x==1,_C);
        e = 0;
        matrix_prior_bal = get_priorbal_system(C,E)
        e = trunc(Int,matrix_prior_bal[i,1]);
        append!(a,e);

        x=0;
        x = rand(findall(x-> ~(x in b) && ~(x in centers_old) && ~(x in centers),CANDIDATAS))
        append!(b,x);

        _C[e] = 0;
        _C[x] = 1;
    end
    return _C;
end

function swap_center_priorbal_grid(C,E,k)
    _C = zeros(Int64,length(CANDIDATAS));
    _C = copy(C);
    centers_old = findall(x->x==1,C)
    a = [];
    b = [];
    for i=1:k
        centers = (x->x==1,_C);
        e = 0;
        matrix_prior_bal = get_priorbal_system(C,E)
        e = trunc(Int,matrix_prior_bal[i,1]);
        append!(a,e);
        zone =get_zone(e);
        x=0
        while true
            x = rand(findall(x->trunc(x)==1,vec(M[zone,:])));
            if (~(x in centers) && ~(x in b) && ~(x in centers_old))
                append!(b,x);
                break;
            end
        end

        _C[e] = 0;
        _C[x] = 1;
    end
    return _C;
end

function connection_calculation()
    num_stations = length(ESTACIONES);
    num_centers = length(CANDIDATAS);
    c = zeros(Int,num_stations,num_centers)

    if instancia == 0
        for i in ESTACIONES
            z1 = get_zone(i);
            for j in CANDIDATAS
                z2 = get_zone(j);
                if (dist[i,j] < dmax) && (same_zone(i,j) || adjacency_matrix[z1,z2] == 1);
                    I_d = 1;
                    c[i,j] = I_d * acc[i,j];
                else
                    c[i,j] = 0;
                end
            end
        end
    else
        for i=1:length(ESTACIONES)
            z1 = get_zone(i);
            for j=1:length(CANDIDATAS)
                z2 = get_zone(j);
                if (dist[i,j] < dmax) && (same_zone(i,j) || adjacency_matrix[z1,z2] == 1);
                    I_d = 1;
                    c[i,j] = I_d * acc[i,j];
                else
                    c[i,j] = 0;
                end
            end
        end
    end
    return c;
end

function validate_connection(_C)
    stations = findall(x->x==0,_C)
    centers  = findall(x->x==1,_C)
    for i=1:length(centers)
        zone = get_zone(centers[i]);
        stations = findall(x->x==1,vec(M[zone,:]));
        sum_aux = 0;
        for j=1:length(stations)
            sum_aux += c[centers[i],stations[j]];
        end
        if sum_aux == 0
            return false;
        end
    end
    return true;
end

function clusterDiameterStats(C, E)
    #distances = zeros(cl,3);
    #distancias = zeros(length(ESTACIONES),1);
    maxDiameterArr = zeros(cl,1);

    iter = 1;
    #iter2 = 1;
    #iter3 = 1;
    for i in findall(x->x==1,C) #Centers loop
        #contStations = 0;
        maxClusterDiameter = 0;
        stations = [];
        #sum = 0;
        stations = get_stations_center(i,E);
        for j in stations
            #sum += dist[i,j];
            #distancias[iter2] = dist[i,j]; #Distancia de cada estacion con su respectivo centro de zona
            #iter2+=1;
            #contStations+=1;
            for s in stations
                if dist[j,s] > maxClusterDiameter
                    maxClusterDiameter = dist[j,s];
                end
            end
        end
        maxDiameterArr[iter] = maxClusterDiameter; #Diametro maximo del cluster de cada zona
        #Diametro de un cluster, maxima distancia entre dos puntos(estaciones) del cluster (de una misma zona)
        #distances[iter,2] = sum; #Distancia sumada de cada estacion con su respectivo centro (separadoen zonas)
        #distances[iter,1] = i;
        #distances[iter,3] = contStations; #Cantidad de estaciones x zona
        iter+=1;
    end

    maxClusterDiam = maximum(maxDiameterArr);
    avgClusterDiam = (sum(maxDiameterArr))/cl;

    #println(maxDiameterArr);
    #println(maxClusterDiam);
    #println(avgClusterDiam);
    #println(distances[:,2]);
    #println(distances[:,1]); 
    #println(distances[:,3]);  

    return maxClusterDiam, avgClusterDiam, maxDiameterArr;
end

function get_prior_system(C,E)
    centers = findall(x->x==1,C)
    pri = [];
    # pri = zeros(cl,3);
    for i=1:length(centers)
        pri = cat(pri,[centers[i] get_prior_center(centers[i],E)]; dims = 1);
    end
    return pri;
end

function get_prior_center(center,E)
    stations_center = get_stations_center(center,E);
    sums_prior = zeros(length(PRIORIDADES));
    prior_array = zeros(length(PRIORIDADES),1)
    for l in PRIORIDADES
        for i in stations_center
            if l == findall(x->x==1,vec(prior[i,:]))[1]
                sums_prior[l] += prior[i,l];
            end
        end
    end

    for i=1:length(sums_prior)
        prior_array[i,1] = abs(floor(length(findall(x->x==1,vec(prior[:,i])))[1]/cl) - sums_prior[i]);
        # prior_array[i,2] = sums_prior[i] - floor(length(findall(x->x==1,prior[:,i]))[1]/cl);
    end
    # result1 = floor(sum(P)/cl) - sum_prior;
    # result2 = sum_prior - floor(sum(P)/cl);
    return prior_array';
end

function get_bal_system(C,E)
    centers = findall(x->x==1,C);
    bal = zeros(cl,2);
    for i=1:length(centers)
        bal[i,1] = centers[i];
        bal[i,2] = abs(get_balance_center(centers[i],E));
    end
    return bal;
end

function get_balance_center(center,E)
    suma_rmas = 0;
    suma_rmenos = 0;
    result = 0;
    stations_center = get_stations_center(center,E)
    for i=1:length(stations_center)
        suma_rmas+=r_mas[stations_center[i]];
        suma_rmenos+=r_menos[stations_center[i]];
    end
    result1 = suma_rmas + suma_rmenos == 0 ? 0 : ((-suma_rmas + suma_rmenos)/(suma_rmas + suma_rmenos));
    result2 = suma_rmas + suma_rmenos == 0 ? 0 : ((suma_rmas - suma_rmenos)/(suma_rmas + suma_rmenos));
    return result1;
    #return result1 <= balance && result2 <= balance;
end

function get_stations_center(center,E)
    x = [];
    if instancia == 0
        for i in ESTACIONES
            if (E[i] == center)
                append!(x,i)
            end
        end
    else
        for i=1:length(ESTACIONES)
            if (E[i] == center)
                append!(x,i)
            end
        end
    end
    return x;
end

function get_distances_cluster_system(C,E)
    distances = zeros(cl,2);
    iter = 1;
    for i in findall(x->x==1,C)
        stations = [];
        sum = 0;
        stations = get_stations_center(i,E);
        for j in stations
            sum += dist[i,j];
        end
        distances[iter,2] = sum;
        distances[iter,1] = i;
        iter+=1;
    end
    return sortslices(distances, dims=1, by=x->(x[2]), rev=true);
end

function get_priorbal_system(C,E)
    a = get_prior_system(C,E);
    r,c = size(a);
    b = get_bal_system(C,E);
    prior_bal_system = zeros(cl,2)
    for i=1:r
        prior_bal_system[i,1] = a[i,1];
        prior_bal_system[i,2] = sum(a[i,2:end]) + b[i,2]*100;
    end
    return prior_bal_system = sortslices(prior_bal_system, dims=1, by=x->(x[2]), rev=true)
end

function get_points(a,b,c,d,M)
    N = size(M)[1];
    points = zeros(length(ESTACIONES));
    index = 1;
    for i=1:N
        if (M[i,3]>a && M[i,3]<=b && M[i,2]<c && M[i,2]>=d)
            # push!(points,trunc(Int64,M[i,1]));
            points[i] = 1;
        end
    end
    return points
end

function get_grid()
    f = open(joinpath(@__DIR__, "geodata.txt"), "r+");
    lines = readlines(f);

    geo_data = zeros(Float64,length(ESTACIONES),3)

    if instancia == 0
        let
            i = 1;
            for line in lines
                data = split(line)
                station = parse(Int64,data[1]);
                latitude = parse(Float64,data[2]);
                longitude = parse(Float64,data[3]);

                geo_data[i,1] = station;
                geo_data[i,2] = latitude;
                geo_data[i,3] = longitude;
                i+=1;
            end
        end
    else
        let
            i = 1;
            for line in lines
                data = split(line)
                station = parse(Int64,data[1]);
                latitude = parse(Float64,data[2]);
                longitude = parse(Float64,data[3]);

                for x=1:length(ESTACIONES)
                    if ESTACIONES[x] == station
                        geo_data[i,1] = station;
                        geo_data[i,2] = latitude;
                        geo_data[i,3] = longitude;
                        i+=1;
                    end
                end
            end
        end
    end
    #lat
    lat_arr = geo_data[:,2];

    #long
    long_arr = geo_data[:,3];


    #min_lat   = minimum(lat_arr);
    min_lat   = 19.3582;
    #max_lat   = maximum(lat_arr);
    max_lat   = 19.444033;

    med_lat  = (max_lat + min_lat) / 2;
    #quarter_lat = (max_lat + med_lat) / 2;
    quarter_lat = 19.42;
    a = min_lat + (med_lat - min_lat)/3;
    b = a + (med_lat - min_lat)/3;

    #min_long  = minimum(long_arr);
    min_long  = -99.20781;
    #max_long  = maximum(long_arr);
    max_long  = -99.13;

    amp1 = (max_long-min_long)/6;
    _1_long = min_long + amp1;
    _2_long = _1_long + amp1;
    _3_long = _2_long + amp1;
    _4_long = _3_long + amp1;
    _5_long = _4_long + amp1;

    c = -99.1910;
    d = -99.15;

    amp2 = (d-c)/3;

    _6_long = c+amp2;
    _7_long = _6_long + amp2;
    amp3 = (d-c)/2;
    _8_long = c + amp3;

    if instancia == 0
        z1 = get_points(min_long,_1_long,max_lat,quarter_lat,geo_data)';
        z2 = get_points(_1_long,_2_long,max_lat,quarter_lat,geo_data)';
        z3 = get_points(_2_long,_3_long,max_lat,quarter_lat,geo_data)';
        z4 = get_points(_3_long,_4_long,max_lat,quarter_lat,geo_data)';
        z5 = get_points(_4_long,_5_long,max_lat,quarter_lat,geo_data)';
        z6 = get_points(_5_long,max_long,max_lat,quarter_lat,geo_data)';
        z7 = get_points(c,_6_long,quarter_lat,med_lat,geo_data)';
        z8 = get_points(_6_long,_7_long,quarter_lat,med_lat,geo_data)';
        z9 = get_points(_7_long,d,quarter_lat,med_lat,geo_data)';
        z10 = get_points(c,_8_long,med_lat,b,geo_data)';
        z11 = get_points(_8_long,d,med_lat,b,geo_data)';
        z12 = get_points(c,_8_long,b,a,geo_data)';
        z13 = get_points(_8_long,d,b,a,geo_data)';
        z14 = get_points(c,_8_long,a,min_lat,geo_data)';
        z15 = get_points(_8_long,d,a,min_lat,geo_data)';

        M = vcat(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15);
        return M;
    elseif instancia == 2
        z1 = get_points(_3_long,_4_long,max_lat,quarter_lat,geo_data)';
        z2 = get_points(_4_long,_5_long,max_lat,quarter_lat,geo_data)'; 
        z3 = get_points(_5_long,max_long,max_lat,quarter_lat,geo_data)'; 
        z4 = get_points(_6_long,d,quarter_lat,med_lat,geo_data)';
        z5 = get_points(_6_long,d,med_lat,b,geo_data)'; 
        z6 = get_points(_6_long,d,b,a,geo_data)'; 
        z7 = get_points(_6_long,d,a,min_lat,geo_data)';

        M = vcat(z1,z2,z3,z4,z5,z6,z7);
        return M;
    else
        z1 = get_points(min_long,_2_long,max_lat,quarter_lat,geo_data)';
        z2 = get_points(_2_long,_4_long,max_lat,quarter_lat,geo_data)';
        z3 = get_points(c,_6_long,quarter_lat,med_lat,geo_data)';
        z4 = get_points(_6_long,_7_long,quarter_lat,med_lat,geo_data)';
        z5 = get_points(c,d,med_lat,b,geo_data)';
        z6 = get_points(c,d,b,a,geo_data)';
        z7 = get_points(c,d,a,min_lat,geo_data)';

        M = vcat(z1,z2,z3,z4,z5,z6,z7);
        return M;
    end 
end

function change_grid()
    min_lat   = 19.3582; #Este valor se mantiene, ya que es el borde inferior de la grilla
    max_lat   = 19.444033; #Este valor se mantiene, ya que es el borde superior de la grilla

    med_lat  = (max_lat + min_lat) / 2; #Este valor se mantiene, ya que es el borde inferior de las zonas 7,8 y 9
    quarter_lat = 19.42; #Este valor se mantiene, ya que es el borde inferior de las primeras 6 zonas
    a = min_lat + (med_lat - min_lat)/3; #Este valor se mantiene, ya que es el borde inferior de las zonas 12 y13
    b = a + (med_lat - min_lat)/3; #Este valor se mantiene, ya que es el borde inferior de las zonas 10 y 11

    min_long  = -99.20781; #Este valor se mantiene, ya que es el borde izquierdo superior de la grilla, 6 primeras zonas
    max_long  = -99.13; #Este valor se mantiene, ya que es el borde derecho superior de la grilla, 6 primeras zonas

    zones_ = []; # ALMACENA ZONAS YA PERTURBADAS ANTERIORMENTE

    z_ = rand(1:15); #Zona aleatoria a perturbar
    t_ = rand(1:2); # 1 = Achicar zona, 2 = Expandir zona
    s_ = rand(1:2); # Sentido de la perturbación de la zona (de frente a la grilla) 1 = Derecha, 2 = Izquierda
    push!(zones_,z_);

    amp1 = (max_long-min_long)/6;

    if z_ == 1
        if t_ == 1
            println("Se achicó la Zona1, se agrandó la Zona2.");
            _1_long = min_long + (amp1/2);
            _2_long = _1_long + (amp1 + amp1/2);
            _3_long = _2_long + amp1;
            _4_long = _3_long + amp1;
            _5_long = _4_long + amp1;
        else
            println("Se agrandó la Zona1, se achicó la Zona2.");
            _1_long = min_long + (amp1 + amp1/2);
            _2_long = _1_long + (amp1/2);
            _3_long = _2_long + amp1;
            _4_long = _3_long + amp1;
            _5_long = _4_long + amp1;
        end
    elseif z_ == 2
        if t_ == 1
            if s_ == 1
                println("Se achicó la Zona2, se agrando la Zona1.");
                _1_long = min_long + (amp1 + amp1/2);
                _2_long = _1_long + (amp1/2);
                _3_long = _2_long + amp1;
                _4_long = _3_long + amp1;
                _5_long = _4_long + amp1;
            else
                println("Se achicó la Zona2, se agrandó la Zona3.");
                _1_long = min_long + amp1;
                _2_long = _1_long + (amp1/2);
                _3_long = _2_long + (amp1 + amp1/2);
                _4_long = _3_long + amp1;
                _5_long = _4_long + amp1;
            end
        else
            if s_ == 1
                println("Se agrandó la Zona2, se achicó la Zona1.");
                _1_long = min_long + (amp1/2);
                _2_long = _1_long + (amp1/2 + amp1);
                _3_long = _2_long + amp1;
                _4_long = _3_long + amp1;
                _5_long = _4_long + amp1;
            else
                println("Se agrandó la Zona2, se achicó la Zona3.");
                _1_long = min_long + amp1;
                _2_long = _1_long + (amp1/2 + amp1);
                _3_long = _2_long + (amp1/2);
                _4_long = _3_long + amp1;
                _5_long = _4_long + amp1;
            end
        end
    elseif z_ == 3
        if t_ == 1
            if s_ == 1
                println("Se achicó la Zona3, se agrandó la Zona2.");
                _1_long = min_long + amp1;
                _2_long = _1_long + (amp1/2 + amp1);
                _3_long = _2_long + (amp1/2);
                _4_long = _3_long + amp1;
                _5_long = _4_long + amp1;
            else
                println("Se achicó la Zona3, se agrandó la Zona4.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + (amp1/2);
                _4_long = _3_long + (amp1/2 + amp1);
                _5_long = _4_long + amp1;
            end
        else
            if s_ == 1
                println("Se agrandó la Zona3, se achicó la Zona2.");
                _1_long = min_long + amp1;
                _2_long = _1_long + (amp1/2);
                _3_long = _2_long + (amp1/2 + amp1);
                _4_long = _3_long + amp1;
                _5_long = _4_long + amp1;
            else
                println("Se agrandó la Zona3, se achicó la Zona4.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + (amp1/2 + amp1);
                _4_long = _3_long + (amp1/2);
                _5_long = _4_long + amp1;
            end
        end
    elseif z_ == 4
        if t_ == 1
            if s_ == 1
                println("Se achicó la Zona4, se agrandó la Zona3.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + (amp1/2 + amp1);
                _4_long = _3_long + (amp1/2);
                _5_long = _4_long + amp1;
            else
                println("Se achicó la Zona4, se agrandó la Zona5.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + amp1;
                _4_long = _3_long + (amp1/2);
                _5_long = _4_long + (amp1/2 + amp1);
            end
        else
            if s_ == 1
                println("Se agrandó la Zona4, se achicó la Zona3.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + (amp1/2);
                _4_long = _3_long + (amp1/2 + amp1);
                _5_long = _4_long + amp1;
            else
                println("Se agrandó la Zona4, se achicó la Zona5.");
                amp1 = (max_long-min_long)/6;
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + amp1;
                _4_long = _3_long + (amp1/2 + amp1);
                _5_long = _4_long + (amp1/2);
            end
        end
    elseif z_ == 5
        if t_ == 1
            if s_ == 1
                println("Se agrandó la Zona5, se achicó la Zona4.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + amp1;
                _4_long = _3_long + (amp1/2);
                _5_long = _4_long + (amp1/2 + amp1);
            else
                println("Se agrandó la Zona5, se achicó la Zona6.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + amp1;
                _4_long = _3_long + amp1;
                _5_long = _4_long + (amp1/2 + amp1);
            end
        else
            if s_ == 1
                println("Se achicó la Zona5, se agrandó la Zona4.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + amp1;
                _4_long = _3_long + (amp1/2 + amp1);
                _5_long = _4_long + (amp1/2);
            else
                println("Se achicó la Zona5, se agrandó la Zona6.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + amp1;
                _4_long = _3_long + amp1;
                _5_long = _4_long + (amp1/2);
            end
        end
    elseif z_ == 6
        if t_ == 1
                println("Se achicó la Zona6, se agrandó la Zona5.");
                _1_long = min_long + amp1;
                _2_long = _1_long + amp1;
                _3_long = _2_long + amp1;
                _4_long = _3_long + amp1;
                _5_long = _4_long + (amp1/2 + amp1);
        else
            println("Se agrandó la Zona6, se achicó la Zona5.");
            _1_long = min_long + amp1;
            _2_long = _1_long + amp1;
            _3_long = _2_long + amp1;
            _4_long = _3_long + amp1;
            _5_long = _4_long + (amp1/2);
        end
    else
        ##SE MANTIENE IGUAL QUE EN LA GRILLA ORIGINAL
        _1_long = min_long + amp1;
        _2_long = _1_long + amp1;
        _3_long = _2_long + amp1;
        _4_long = _3_long + amp1;
        _5_long = _4_long + amp1;
    end

    c = -99.1910; #Este valor se mantiene, ya que es el borde izquierdo inferior de la grilla, 9 siguientes zonas
    d = -99.15; #Este valor se mantiene, ya que es el borde derecho inferior de la grilla, 9 siguientes zonas

    amp2 = (d-c)/3;

    if z_ == 7
        if t_ == 1
            println("Se achicó la Zona7, se agrandó la Zona8.");
            _6_long = c+(amp2/2);
            _7_long = _6_long + (amp2/2 + amp2);
        else
            println("Se agrandó la Zona7, se achicó la Zona8.");
            _6_long = c+(amp2/2 + amp2);
            _7_long = _6_long + (amp2/2);
        end
    elseif z_ == 8
        if t_ == 1
            if s_ == 1
                println("Se achicó la Zona8, se agrandó la Zona7.");
                _6_long = c+(amp2/2 + amp2);
                _7_long = _6_long + (amp2/2);
            else
                println("Se achicó la Zona8, se agrandó la Zona9.");
                _6_long = c+(amp2);
                _7_long = _6_long + (amp2/2);
            end
        else
            if s_ == 1
                println("Se agrandó la Zona8, se achicó la Zona7.");
                _6_long = c+(amp2/2);
                _7_long = _6_long + (amp2/2 + amp2);
            else
                println("Se agrandó la Zona8, se achicó la Zona9.");
                _6_long = c+(amp2);
                _7_long = _6_long + (amp2/2 + amp2);
            end
        end
    elseif z_ == 9
        if t_ == 1
            println("Se achicó la Zona9, se agrandó la Zona8.");
            _6_long = c+amp2;
            _7_long = _6_long + (amp2/2 + amp2);
        else
            println("Se agrandó la Zona9, se achicó la Zona8.");
            _6_long = c+amp2;
            _7_long = _6_long + (amp2/2);
        end
    else
        ##SE MANTIENEN LOS VALORES ORIGINALES
        _6_long = c+amp2;
        _7_long = _6_long + amp2;
    end

    z1 = get_points(min_long,_1_long,max_lat,quarter_lat,geo_data)';
    z2 = get_points(_1_long,_2_long,max_lat,quarter_lat,geo_data)';
    z3 = get_points(_2_long,_3_long,max_lat,quarter_lat,geo_data)';
    z4 = get_points(_3_long,_4_long,max_lat,quarter_lat,geo_data)';
    z5 = get_points(_4_long,_5_long,max_lat,quarter_lat,geo_data)';
    z6 = get_points(_5_long,max_long,max_lat,quarter_lat,geo_data)';
    z7 = get_points(c,_6_long,quarter_lat,med_lat,geo_data)';
    z8 = get_points(_6_long,_7_long,quarter_lat,med_lat,geo_data)';
    z9 = get_points(_7_long,d,quarter_lat,med_lat,geo_data)';

    amp3 = (d-c)/2;
    _8_long = c + amp3;

    if z_ == 10
        if t_ == 1
            println("Se achicó la Zona10, se agrandó la Zona11.");
            z10 = get_points(c,(_8_long - (amp3/2)),med_lat,b,geo_data)';
            z11 = get_points((_8_long + (amp3/2)),d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        else
            println("Se agrandó la Zona10, se achicó la Zona11.");
            z10 = get_points(c,(_8_long + (amp3/2)),med_lat,b,geo_data)';
            z11 = get_points((_8_long - (amp3/2)),d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        end
    elseif z_ == 11
        if t_ == 1
            println("Se achicó la Zona11, se agrandó la Zona10.");
            z10 = get_points(c,(_8_long + (amp3/2)),med_lat,b,geo_data)';
            z11 = get_points((_8_long - (amp3/2)),d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        else
            println("Se agrandó la Zona11, se achicó la Zona10.");
            z10 = get_points(c,(_8_long - (amp3/2)),med_lat,b,geo_data)';
            z11 = get_points((_8_long + (amp3/2)),d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        end
    elseif z_ == 12
        if t_ == 1
            println("Se achicó la Zona12, se agrandó la Zona13.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,(_8_long - (amp3/2)),b,a,geo_data)';
            z13 = get_points((_8_long + (amp3/2)),d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        else
            println("Se agrandó la Zona12, se achicó la Zona13.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,(_8_long + (amp3/2)),b,a,geo_data)';
            z13 = get_points((_8_long - (amp3/2)),d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        end
    elseif z_ == 13
        if t_ == 1
            println("Se achicó la Zona13, se agrandó la Zona12.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,(_8_long + (amp3/2)),b,a,geo_data)';
            z13 = get_points((_8_long - (amp3/2)),d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        else
            println("Se agrandó la Zona13, se achicó la Zona12.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,(_8_long - (amp3/2)),b,a,geo_data)';
            z13 = get_points((_8_long + (amp3/2)),d,b,a,geo_data)';
            z14 = get_points(c,_8_long,a,min_lat,geo_data)';
            z15 = get_points(_8_long,d,a,min_lat,geo_data)';
        end
    elseif z_ == 14
        if t_ == 1
            println("Se achicó la Zona14, se agrandó la Zona15.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,(_8_long - (amp3/2)),a,min_lat,geo_data)';
            z15 = get_points((_8_long + (amp3/2)),d,a,min_lat,geo_data)';
        else
            println("Se agrandó la Zona14, se achicó la Zona15.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,(_8_long + (amp3/2)),a,min_lat,geo_data)';
            z15 = get_points((_8_long - (amp3/2)),d,a,min_lat,geo_data)';
        end
    elseif z_ == 15
        if t_ == 1
            println("Se achicó la Zona15, se agrandó la Zona14.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,(_8_long + (amp3/2)),a,min_lat,geo_data)';
            z15 = get_points((_8_long - (amp3/2)),d,a,min_lat,geo_data)';
        else
            println("Se agrandó la Zona15, se achicó la Zona14.");
            z10 = get_points(c,_8_long,med_lat,b,geo_data)';
            z11 = get_points(_8_long,d,med_lat,b,geo_data)';
            z12 = get_points(c,_8_long,b,a,geo_data)';
            z13 = get_points(_8_long,d,b,a,geo_data)';
            z14 = get_points(c,(_8_long - (amp3/2)),a,min_lat,geo_data)';
            z15 = get_points((_8_long + (amp3/2)),d,a,min_lat,geo_data)';
        end
    else
        z10 = get_points(c,_8_long,med_lat,b,geo_data)';
        z11 = get_points(_8_long,d,med_lat,b,geo_data)';
        z12 = get_points(c,_8_long,b,a,geo_data)';
        z13 = get_points(_8_long,d,b,a,geo_data)';
        z14 = get_points(c,_8_long,a,min_lat,geo_data)';
        z15 = get_points(_8_long,d,a,min_lat,geo_data)';
    end

    #NUEVA GRILLA, CON UNA ZONA PERTURBADA
    M = vcat(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15);
    #return M;
end

function same_zone(x,y)
    flag =false;
    zone1 = get_zone(x);
    zone2 = get_zone(y);
    if zone1 == zone2
        flag = true;
    end

    return flag;
end

function get_zone(x)
    return findall(x->trunc(x)==1,M[:,x])[1];
end

function get_adjacency_matrix()
    if instancia == 0
        adjacency_matrix = zeros(Int,15,15)
        adjacency_matrix[1,1] = 1;
        adjacency_matrix[1,2] = 1;
        adjacency_matrix[2,1] = 1;
        adjacency_matrix[2,2] = 1;
        adjacency_matrix[2,3] = 1;
        adjacency_matrix[2,7] = 1;
        adjacency_matrix[3,2] = 1;
        adjacency_matrix[3,3] = 1;
        adjacency_matrix[3,4] = 1;
        adjacency_matrix[3,7] = 1;
        adjacency_matrix[3,8] = 1;
        adjacency_matrix[4,3] = 1;
        adjacency_matrix[4,4] = 1;
        adjacency_matrix[4,5] = 1;
        #adjacency_matrix[4,7] = 1; ##PREGUNTAR ESTO
        adjacency_matrix[4,8] = 1;
        adjacency_matrix[4,9] = 1;
        adjacency_matrix[5,4] = 1;
        adjacency_matrix[5,5] = 1;
        adjacency_matrix[5,6] = 1;
        adjacency_matrix[5,9] = 1;
        adjacency_matrix[6,5] = 1;
        adjacency_matrix[6,6] = 1;
        #adjacency_matrix[6,9] = 1; ##PREGUNTAR ESTO
        adjacency_matrix[7,2] = 1;
        adjacency_matrix[7,3] = 1;
        #adjacency_matrix[7,4] = 1; ##PREGUNTAR ESTO
        adjacency_matrix[7,7] = 1;
        adjacency_matrix[7,8] = 1;
        adjacency_matrix[7,10] = 1;
        adjacency_matrix[8,3] = 1;
        adjacency_matrix[8,4] = 1;
        adjacency_matrix[8,7] = 1;
        adjacency_matrix[8,8] = 1;
        adjacency_matrix[8,9] = 1;
        adjacency_matrix[8,10] = 1;
        adjacency_matrix[8,11] = 1;
        adjacency_matrix[9,4] = 1;
        adjacency_matrix[9,5] = 1;
        #adjacency_matrix[9,6] = 1; ##PREGUNTAR ESTO
        adjacency_matrix[9,8] = 1;
        adjacency_matrix[9,9] = 1;
        adjacency_matrix[9,11] = 1;
        adjacency_matrix[10,7] = 1;
        adjacency_matrix[10,8] = 1;
        adjacency_matrix[10,10] = 1;
        adjacency_matrix[10,11] = 1;
        adjacency_matrix[10,12] = 1;
        adjacency_matrix[10,13] = 1;
        adjacency_matrix[11,8] = 1;
        adjacency_matrix[11,9] = 1;
        adjacency_matrix[11,10] = 1;
        adjacency_matrix[11,11] = 1;
        adjacency_matrix[11,12] = 1;
        adjacency_matrix[11,13] = 1;
        adjacency_matrix[12,10] = 1;
        adjacency_matrix[12,11] = 1;
        adjacency_matrix[12,12] = 1;
        adjacency_matrix[12,13] = 1;
        adjacency_matrix[12,14] = 1;
        adjacency_matrix[12,15] = 1;
        adjacency_matrix[13,10] = 1;
        adjacency_matrix[13,11] = 1;
        adjacency_matrix[13,12] = 1;
        adjacency_matrix[13,13] = 1;
        adjacency_matrix[13,14] = 1;
        adjacency_matrix[13,15] = 1;
        adjacency_matrix[14,12] = 1;
        adjacency_matrix[14,13] = 1;
        adjacency_matrix[14,14] = 1;
        adjacency_matrix[14,15] = 1;
        adjacency_matrix[15,12] = 1;
        adjacency_matrix[15,13] = 1;
        adjacency_matrix[15,14] = 1;
        adjacency_matrix[15,15] = 1;
        return adjacency_matrix;
    elseif instancia == 2
        adjacency_matrix = zeros(Int,15,15)
        adjacency_matrix[1,1] = 1;
        adjacency_matrix[1,2] = 1;
        adjacency_matrix[1,4] = 1;
        adjacency_matrix[2,1] = 1;
        adjacency_matrix[2,2] = 1;
        adjacency_matrix[2,3] = 1;
        adjacency_matrix[2,4] = 1;
        adjacency_matrix[3,2] = 1;
        adjacency_matrix[3,3] = 1;
        adjacency_matrix[4,1] = 1;
        adjacency_matrix[4,2] = 1;
        adjacency_matrix[4,4] = 1;
        adjacency_matrix[4,5] = 1;
        adjacency_matrix[5,4] = 1;
        adjacency_matrix[5,5] = 1;
        adjacency_matrix[5,6] = 1;
        adjacency_matrix[6,5] = 1;
        adjacency_matrix[6,6] = 1;
        adjacency_matrix[6,7] = 1;
        adjacency_matrix[7,6] = 1;
        adjacency_matrix[7,7] = 1;
        return adjacency_matrix;
    else
        adjacency_matrix = zeros(Int,15,15)
        adjacency_matrix[1,1] = 1;
        adjacency_matrix[1,2] = 1;
        adjacency_matrix[1,3] = 1;
        adjacency_matrix[2,1] = 1;
        adjacency_matrix[2,2] = 1;
        adjacency_matrix[2,3] = 1;
        adjacency_matrix[2,4] = 1;
        adjacency_matrix[3,1] = 1;
        adjacency_matrix[3,2] = 1;
        adjacency_matrix[3,3] = 1;
        adjacency_matrix[3,4] = 1;
        adjacency_matrix[3,5] = 1;
        adjacency_matrix[4,2] = 1;
        adjacency_matrix[4,3] = 1;
        adjacency_matrix[4,4] = 1;
        adjacency_matrix[4,5] = 1;
        adjacency_matrix[5,3] = 1;
        adjacency_matrix[5,4] = 1;
        adjacency_matrix[5,5] = 1;
        adjacency_matrix[5,6] = 1;
        adjacency_matrix[6,5] = 1;
        adjacency_matrix[6,6] = 1;
        adjacency_matrix[6,7] = 1;
        adjacency_matrix[7,6] = 1;
        adjacency_matrix[7,7] = 1;
        return adjacency_matrix;
    end
end
