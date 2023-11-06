function [S, V] = compute_geometry_S_V(file_path)
    % WE SHOULD CHECK IF SEGMENT CONTAINS SEGMENT; NOT IF EQUAL

    % Read the data from the text file
    data = dlmread(file_path);
    
    % Extract the columns for position and dimensions
    x_pos = data(:, 1);
    y_pos = data(:, 2);
    z_pos = data(:, 3);
    x_length = data(:, 4);
    y_length = data(:, 5);
    z_length = data(:, 6);
    
    % Initialize total surface area and volume
    total_perimeter = 0;
    total_surface_area = 0;

    % Create a cell array to store the segments
    segments = {};
    
    % Create an array to store the results of check_segment_equal
    segmentEqualResults = {};

    % Loop through each rectangle
    for i = 1:length(x_pos)
        % Calculate the surface area of the current rectangle
        surface_area = x_length(i) * y_length(i);
        
        % Add the current surface area to the total
        total_surface_area = total_surface_area + surface_area;

        % Calculate the perimeter of the current rectangle
        perimeter = 2 * (x_length(i) + y_length(i));

        % Add the current perimeter to the total
        total_perimeter = total_perimeter + perimeter;


        % Define and add the coordinates of the individual segments one by one
        top_segment = [x_pos(i), y_pos(i); x_pos(i) + x_length(i), y_pos(i)]; % Top
        
        segmentEqualResult = check_segment_equal(segments, top_segment);
        segmentEqualResults{end + 1} = segmentEqualResult;
        segments{end + 1} = top_segment;

        if segmentEqualResult == true
            total_perimeter = total_perimeter - 2*x_length(i);
        end

        plot(top_segment(:,1), top_segment(:,2));
        
        
        right_segment = [x_pos(i) + x_length(i), y_pos(i); x_pos(i) + x_length(i), y_pos(i) + y_length(i)]; % Right
        
        segmentEqualResult = check_segment_equal(segments, right_segment);
        segmentEqualResults{end + 1} = segmentEqualResult;
        segments{end + 1} = right_segment;

        if segmentEqualResult == true
            total_perimeter = total_perimeter - 2*y_length(i);
        end

        plot(right_segment(:,1), right_segment(:,2));
        

        bottom_segment = [x_pos(i) + x_length(i), y_pos(i) + y_length(i); x_pos(i), y_pos(i) + y_length(i)]; % Bottom
        
        segmentEqualResult = check_segment_equal(segments, bottom_segment);
        segmentEqualResults{end + 1} = segmentEqualResult;
        segments{end + 1} = bottom_segment;
        
        if segmentEqualResult == true
            total_perimeter = total_perimeter - 2*x_length(i);
        end

        plot(bottom_segment(:,1), bottom_segment(:,2));
        

        left_segment = [x_pos(i), y_pos(i) + y_length(i); x_pos(i), y_pos(i)]; % Left
        
        segmentEqualResult = check_segment_equal(segments, left_segment);
        segmentEqualResults{end + 1} = segmentEqualResult;
        segments{end + 1} = left_segment;

        if segmentEqualResult == true
            total_perimeter = total_perimeter - 2*y_length(i);
        end

        plot(left_segment(:,1), left_segment(:,2));
        

        hold on;
    end
    
    % Set axis labels and title for the 2D plot
    xlabel('X');
    ylabel('Y');
    title('Rectangle Segments');
    
    % Return the results
    S = total_perimeter * z_length(1) + 2 * total_surface_area;
    V = total_surface_area * z_length(1);
end