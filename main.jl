global instancia = args[9]; #0 = PT, 1 = PT_V1, 2 = PT_V2
global last_obj = 0;

include("load_data.jl");
include("gurobi.jl");
include("helpers.jl");
include("RVNS.jl");

using Statistics;

#Variables globales de balance y prioridad.
#Restricción de all prior # allprior = 0, allprior = 1;
#Criterio  # 1 = swap center random,  2 = swap center max distance,  3 = swap center prior bal

#Grilla
global M                = get_grid();

#Matriz de adyacencia de zonas
global adjacency_matrix = get_adjacency_matrix();

#Matriz de conexiones
global c                = connection_calculation();

#parámetros
balance   = args[3];
prioridad = args[4];
dmax      = args[7];
allprior  = args[8];

#= MAIN =#
# Variables metaheurística
r_max = args[2]; #Número iteraciones.
len_N = args[6];  #Tamaño de los vecindarios.
neighborhood_structure = args[5]; #Tamaños estructuras de entorno.
k_max = length(neighborhood_structure); #k máximo.

#Numero de experimentos a realizar.
experimentos = args[1];

println("Parámetros");
println("   experimentos   ",args[1]);
println("   iteraciones    ",args[2]);
println("   balance        ",args[3]);
println("   prioridad      ",args[4]);
println("   mov_vecindario ",args[5]);
println("   tam_vecindario ",args[6]);
println("   dmax           ",args[7]);
println(" ");

#Limite de no mejoras.
const NO_IMPROVE_LIMIT = (r_max / 4);

objs_iter = 0;
objs_array = [];
exp_time_array = [];
if criterio == 1
    println("RVNS con Random Swap Center.")
elseif criterio == 2
    println("RVNS con Max Distance Swap Center.")
else
    println("RVNS con Prior Bal Swap Center.")
end
for e = 1:experimentos
    C_test = zeros(Int64,length(CANDIDATAS));
    E_test = zeros(Int64,length(ESTACIONES));
    exp_time = @elapsed C_test,E_test,objs_iter,improve,obj = @time RVNS(k_max,r_max,len_N,neighborhood_structure,e,NO_IMPROVE_LIMIT,criterio);
    append!(objs_array,objs_iter);
    append!(exp_time_array, exp_time);
    name = "$(balance)_$(prioridad)_exp_$(e)_$(r_max)_$(len_N)_$(improve)_$(obj)";
    filename = name*".txt" 
    open(joinpath("C:/Users/joaqu/Desktop/RVNS final/Resultados finales/$(nombre_instancia)/$(balance)_$(prioridad)/$(r_max)/$(neighborhood_structure)_$(len_N)", filename), "a") do file
        write(file, "tiempo       = $exp_time \n");
    end
end
let suma = 0.0, sumatimes = 0.0;
    for i=1:length(objs_array)
        suma = suma + objs_array[i];
    end
    for i=1:length(exp_time_array)
        sumatimes = sumatimes + exp_time_array[i];
    end

    time_prom = sumatimes/length(exp_time_array);
    promedio  = suma/length(objs_array);
    de    = std(floor.(objs_array));
    best  = minimum(objs_array);
    worst = maximum(objs_array);

    #Resumen resultados
    name = "result_exp_$(balance)_$(prioridad)_$(experimentos)_$(best)";
    filename = name*".txt"
    open(joinpath("C:/Users/joaqu/Desktop/RVNS final/Resultados finales/$(nombre_instancia)/$(balance)_$(prioridad)/$(r_max)/$(neighborhood_structure)_$(len_N)", filename), "w") do file
        write(file, "experimentos   = $experimentos \n");
        write(file, "promedio       = $promedio \n");
        write(file, "d.e            = $de   \n");
        write(file, "best           = $best \n");
        write(file, "worst          = $worst \n");
        write(file, "tiempo_prom    = $time_prom s\n");
        write(file, "tiempo_c/exp   = $exp_time_array\n");
        if criterio == 1
            write(file, "criterio       = swap center random \n");
        elseif criterio == 2
            write(file, "criterio       = swap center max distance \n");
        else
            write(file, "criterio       = swap center prior bal \n");
        end
    end
end
