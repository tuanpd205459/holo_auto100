function zernike_reconstruction(coefficients, real_radius, lambda, scale_z)
    % Function to reconstruct a surface from Zernike coefficients
    % coefficients: Vector of Zernike coefficients
    % real_radius: Radius of the circular aperture in mm
    % lambda: Wavelength in nm
    % scale_z: Scaling factor for Z-axis (e.g., 0.2 to shrink Z)

    % Number of coefficients
    num_coeff = length(coefficients);

    % Generate a grid of points within the actual radius
    [x, y] = meshgrid(linspace(-real_radius, real_radius, 400), ...
                       linspace(-real_radius, real_radius, 400));
    [theta, rho] = cart2pol(x, y);
    
    % Normalize rho to [0, 1] for Zernike polynomial calculation
    rho = rho / real_radius;

    % Mask for points within the actual circle
    mask = rho <= 1;

    % Initialize the surface
    surface = zeros(size(x));

    % Calculate the surface using the Zernike coefficients
    n = 0; % Radial degree
    m = 0; % Azimuthal frequency
    for k = 1:num_coeff
        % Calculate the Zernike polynomial
        Z = zernike_polynomial(n, m, rho(mask), theta(mask));
        
        % Add the contribution of the current Zernike polynomial
        surface(mask) = surface(mask) + coefficients(k) * Z;

        % Update the radial degree and azimuthal frequency
        if m == n
            n = n + 1;
            m = -n;
        else
            m = m + 2;
        end
    end

    % Scale surface height (reduce Z exaggeration)
    surface = (surface / lambda) * scale_z; % Shrink Z-axis

    % Set values outside the mask to NaN
    surface(~mask) = NaN;

    % Plot the reconstructed surface
    figure;
    surf(x, y, surface, 'EdgeColor', 'none');

    % Formatting for better visualization
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel(['Surface Height (', num2str(scale_z), 'λ)']); % Z hiển thị có scale
    title(['Reconstructed Surface (λ = ', num2str(lambda), ' nm)']);
    axis equal; % Giữ tỷ lệ đúng giữa các trục X, Y
    daspect([1 1 scale_z]); % Điều chỉnh tỷ lệ trực quan X:Y:Z
    colormap jet;
    colorbar;
    caxis([-1 1] * max(abs(surface(:)))); % Cân bằng thang màu
    view(3);
end



function Z = zernike_polynomial(n, m, rho, theta)
    % Function to calculate the Zernike polynomial Z_n^m(rho, theta)
    R = zernike_radial(n, abs(m), rho);
    if m >= 0
        Z = R .* cos(m * theta);
    else
        Z = R .* sin(abs(m) * theta);
    end
end

function R = zernike_radial(n, m, rho)
    % Function to calculate the radial part of the Zernike polynomial
    R = zeros(size(rho));
    for s = 0:((n - m) / 2)
        c = (-1)^s * factorial(n - s) / (factorial(s) * factorial((n + m) / 2 - s) * factorial((n - m) / 2 - s));
        R = R + c * rho.^(n - 2 * s);
    end
end
