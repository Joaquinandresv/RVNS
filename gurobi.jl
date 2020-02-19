using JuMP, Gurobi

#Contador infactibles
global cont_infactible  = 0;

function Gurobi_optimal(C)
    num_stations  = length(ESTACIONES);
    num_candidatas  = length(CANDIDATAS);
    E = zeros(Int64,num_stations);
    m = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0)) #No muestra resultados por consola.

    @variable(m,x[i=1:num_stations,j=1:num_candidatas],Bin)

    if instancia == 0
        @objective(m,Min,sum(dist[i, j] * x[i, j] for i in ESTACIONES, j in CANDIDATAS));
        for i in ESTACIONES
            @constraint(m,sum(x[i,j] for j in CANDIDATAS) == 1)
        end

        for i in ESTACIONES
            for j in CANDIDATAS
                @constraint(m,x[i,j] <= c[i,j]*C[j])
            end
        end

        for j in CANDIDATAS 
            if C[j] == 1
                if allprior == 1
                    for l in PRIORIDADES
                        pxsum = @expression(m, sum(prior[i,l]*x[i,j] for i in ESTACIONES))
                        psum = @expression(m, sum(prior[i,l] for i in ESTACIONES))
                        @constraint(m,(pxsum  - floor(psum/cl)*C[j]) <= prioridad)
                        @constraint(m,(floor(psum/cl)*C[j] - pxsum) <= prioridad)
                    end 
                else
                    l = 1
                    pxsum = @expression(m, sum(prior[i,l]*x[i,j] for i in ESTACIONES))
                    psum = @expression(m, sum(prior[i,l] for i in ESTACIONES))
                    @constraint(m,(pxsum  - floor(psum/cl)*C[j]) <= prioridad)
                    @constraint(m,(floor(psum/cl)*C[j] - pxsum) <= prioridad)
                end
            end
        end

        for j in CANDIDATAS
            expr2 = @expression(m,sum(r_menos[i]*x[i,j] for i in ESTACIONES));
            expr3 = @expression(m,sum(r_mas[i]*x[i,j] for i in ESTACIONES));
            @constraint(m, (expr2 - expr3)  <= balance  *  (expr2 + expr3));
            @constraint(m,(-expr2 + expr3) <= balance  *  (expr2 + expr3));
        end

        optimize!(m);
        status = termination_status(m);
        if status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("INFACTIBLE, BUSCANDO OTRA SOLUCIÓN");
            #cont_infactible = cont_infactible + 1;
            return Inf, zeros(num_stations,num_stations); 
        end

        #println("TERMINATION STATUS: ", status);
        Z_opt = objective_value(m);
        #println("FUNCION OBJETIVO POR RETORNAR: ", round(Z_opt));
        x_opt = value.(x);
        #println("X_OPT: ", length(x_opt));

        if (status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED) || (round(Z_opt) - floor(round(Z_opt)) != 0 || length(x_opt) == 0)
            return Inf, zeros(num_stations,num_stations);
        else
            for i in ESTACIONES
                arrayE = findall(x->x==1,x_opt[i,:]);
                if length(arrayE) > 0
                    E[i] = arrayE[1];
                end
            end
            return round(Z_opt), E
        end
    else
        @objective(m,Min,sum(dist[i, j] * x[i, j] for i=1:length(ESTACIONES), j=1:length(CANDIDATAS)));
        for i=1:length(ESTACIONES)
            @constraint(m,sum(x[i,j] for j=1:length(CANDIDATAS)) == 1)
        end

        for i=1:length(ESTACIONES)
            for j=1:length(CANDIDATAS)
                @constraint(m,x[i,j] <= c[i,j]*C[j])
            end
        end

        for j=1:length(CANDIDATAS) 
            if C[j] == 1
                if allprior == 1
                    for l in PRIORIDADES
                        pxsum = @expression(m, sum(prior[i,l]*x[i,j] for i=1:length(ESTACIONES)))
                        psum = @expression(m, sum(prior[i,l] for i=1:length(ESTACIONES)))
                        @constraint(m,(pxsum  - floor(psum/cl)*C[j]) <= prioridad)
                        @constraint(m,(floor(psum/cl)*C[j] - pxsum) <= prioridad)
                    end 
                else
                    l = 1
                    pxsum = @expression(m, sum(prior[i,l]*x[i,j] for i=1:length(ESTACIONES)))
                    psum = @expression(m, sum(prior[i,l] for i=1:length(ESTACIONES)))
                    @constraint(m,(pxsum  - floor(psum/cl)*C[j]) <= prioridad)
                    @constraint(m,(floor(psum/cl)*C[j] - pxsum) <= prioridad)
                end
            end
        end

        for j=1:length(CANDIDATAS)
            expr2 = @expression(m,sum(r_menos[i]*x[i,j] for i=1:length(ESTACIONES)));
            expr3 = @expression(m,sum(r_mas[i]*x[i,j] for i=1:length(ESTACIONES)));
            @constraint(m, (expr2 - expr3)  <= balance  *  (expr2 + expr3));
            @constraint(m,(-expr2 + expr3) <= balance  *  (expr2 + expr3));
        end

        optimize!(m);
        status = termination_status(m);
        if status == MOI.INFEASIBLE_OR_UNBOUNDED
            #ultimo solución objetivo
            last_obje        = 0;
            println("INFACTIBLE, BUSCANDO OTRA SOLUCIÓN");
            if cont_infactible == 0
                global cont_infactible = cont_infactible + 1;
                return Inf, zeros(num_stations,num_stations); 
            elseif (cont_infactible == 0 || cont_infactible > 0) && last_obje != 0
                global cont_infactible = cont_infactible + 1;
                return last_obje, zeros(num_stations,num_stations);
            else
                global cont_infactible = cont_infactible + 1;                
                return Inf, zeros(num_stations,num_stations);
            end
        end

        #println("TERMINATION STATUS: ", status);
        Z_opt = objective_value(m);
        #println("FUNCION OBJETIVO POR RETORNAR: ", round(Z_opt));
        x_opt = value.(x);
        global last_obj = round(Z_opt);
        last_obje = last_obj;
        #println("X_OPT: ", length(x_opt));

        if (status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED) || (round(Z_opt) - floor(round(Z_opt)) != 0 || length(x_opt) == 0)
            return Inf, zeros(num_stations,num_stations);
        else
            for i=1:length(ESTACIONES)
                arrayE = findall(x->x==1,x_opt[i,:]);
                if length(arrayE) > 0
                    E[i] = arrayE[1];
                end
            end
            return round(Z_opt), E
        end    
    end
end
