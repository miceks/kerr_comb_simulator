function batch = sim_seq(batch, data_queue)
% Simulate batch using sequential for loop

    % Reorder systems to make active ones contiguous
    [~, order_idx] = sort([batch.systems.budget] .* ~[batch.systems.frozen], 'descend');
    systems = batch.systems(order_idx);
    last = find([systems.budget] > 0 & ~[systems.frozen], 1, 'last');
    if isempty(last)
        last = 0;
    end


    for idx = 1:last

        % Exhaust budget of system
        for series = 1:systems(idx).budget
            
            switch batch.type
                case 'Ikeda', systems(idx) = split_step_Ikeda(systems(idx));
                case 'LLE',   systems(idx) = split_step_LLE(systems(idx));
            end

            % Update state and freeze counter
            systems(idx) = update_state(systems(idx), systems(idx).series, batch.conv_cw_tol, batch.conv_var_tol);

            if systems(idx).state(end) == systems(idx).state(end - 1)
                systems(idx).freeze_in = systems(idx).freeze_in - 1;
            else
                systems(idx).freeze_in = systems(idx).freeze_lim;
            end
        
            if systems(idx).freeze_in <= 0
                systems(idx).frozen = true;
                break
            end

        end

        % Broadcast progress
        send(data_queue, systems(idx).budget);
        systems(idx).budget = 0;

    end

    % Reset previous order of systems
    [~, order_idx] = sort(order_idx, 'ascend');
    batch.systems = systems(order_idx);

end