function addReserve(model, coord, optimal_force, max_control)
    
    %Short function to add reserves to model
    % Inputs:
        % model - OpenSim model object
        % coord - string of coordinate to add reserve to
        % optimal force - optimal force value for actuator
        % max_control - maximum control for actuator

    %Import opensim libraries
    import org.opensim.modeling.* 
       
    actu = ActivationCoordinateActuator();
    actu.set_coordinate(coord);
    if startsWith(coord,'pelvis')
        prefix = 'residual_';
    else
        prefix = 'reserve_';
    end
    actu.setName([prefix,coord]);
    actu.setOptimalForce(optimal_force);
    actu.setMinControl(-max_control);
    actu.setMaxControl(max_control);
    model.addForce(actu);
    
end