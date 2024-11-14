% Paramètres pour l'extension à 15%
ex1 = 0.08;
ey1 = 0.024;
exy1 = 0.024;

% Calcul des déformations normales et de cisaillement pour l'extension à 15%
Phi = linspace(0, 360, 180); % Angles de 0° à 360°
PhiRad = Phi * pi / 180; % Conversion en radians
en1 = (ex1 + ey1) / 2 + ((ex1 - ey1) / 2) * cos(2 * PhiRad) + exy1 * sin(2 * PhiRad);
yn1 = (ex1 - ey1) / 2 * sin(2 * PhiRad) + exy1 * cos(2 * PhiRad);

% Paramètres pour l'extension à 10%
ex2 = 0.146;
ey2 = 0.042;
exy2 = 1.48;

% Calcul des déformations normales et de cisaillement pour l'extension à 10%
en2 = (ex2 + ey2) / 2 + ((ex2 - ey2) / 2) * cos(2 * PhiRad) + exy2 * sin(2 * PhiRad);
yn2 = (ex2 - ey2) / 2 * sin(2 * PhiRad) + exy2 * cos(2 * PhiRad);

% Tracé des deux extensions sur la même figure avec des sous-graphiques
figure;

% Tracé pour l'extension à 15%
subplot(1, 2, 1); % Première position
polarplot(PhiRad, en1, 'b', 'LineWidth', 2);
hold on;
polarplot(PhiRad, yn1, 'r--', 'LineWidth', 2);
title('Extension à 15%');
legend('\epsilon (Déformation normale)', '\gamma (Distorsion)');

% Tracé pour l'extension à 10%
subplot(1, 2, 2); % Deuxième position
polarplot(PhiRad, en2, 'b', 'LineWidth', 2);
hold on;
polarplot(PhiRad, yn2, 'r--', 'LineWidth', 2);
title('Extension à 10%');
legend('\epsilon (Déformation normale)', '\gamma (Distorsion)');
