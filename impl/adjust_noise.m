function batch = adjust_noise(batch, id_range, noise)
% Change noise of systems

    if all(id_range == 0)
        id_range = 1:length(batch.systems);
    end

    for id = id_range
        batch.systems(id).noise(end) = noise;
    end

end