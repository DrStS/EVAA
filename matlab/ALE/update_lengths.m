function [lengths] = update_lengths(original_spring_length, position_vector, corners, angles)
    C = get_corners(corners, angles);
    lengths(1) = original_spring_length(1) + C(3,1) + position_vector(1) - position_vector(4);
    lengths(2) = original_spring_length(2) + position_vector(4) - position_vector(5);
    lengths(3) = original_spring_length(3) + C(3,2) + position_vector(1) - position_vector(6);
    lengths(4) = original_spring_length(4) + position_vector(6) - position_vector(7);
    lengths(5) = original_spring_length(5) + C(3,3) + position_vector(1) - position_vector(8);
    lengths(6) = original_spring_length(6) + position_vector(8) - position_vector(9);
    lengths(7) = original_spring_length(7) + C(3,4) + position_vector(1) - position_vector(10);
    lengths(8) = original_spring_length(8) + position_vector(10) - position_vector(11);
end