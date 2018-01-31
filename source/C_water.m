 function rhocp = C_water(T)

% Finds the specific heat and volumetric heat capacity of liquid water 
% within the temperature range, 236-360 K.

% Notation:

%   T     = temperature (C)                (vector)
%   rhocp = volumetric heat capacity       (vector)

% Options:

%   dopt  = 1 display a flag if Tk < 236 K and pause
%         = 2 use the C_water(236 K) value at colder temps and keep running

 dopt = 2;
% ______________________________________________

% Written by:

%   Gary Clow
%   Institute of Arctic and Alpine Research
%   University of Colorado
%   Boulder, Colorado USA
%   Email: gary.clow@colorado.edu
% ______________________________________________

 Tk = T + 273.15;   % kelvin temperature

% parameters

 rhow = 1000;       % density of water 

% table of Cp_water(Tk) values

 A = [236   5699; 237 5501.2; 238 5351.3; 239 5228.7; 240 5123.5; ...
      241 5030.8; 242 4947.8; 243 4872.8; 244 4805.1; 245 4743.7; ...
      246 4688.1; 247 4637.8; 248 4592.4; 249 4551.3; 250 4514.3; ...
      251   4481; 252   4451; 253   4424; 254 4399.8; 255   4378; ...
      256 4358.4; 257 4340.8; 258 4325.1; 259 4310.9; 260 4298.2; ...
      261 4286.8; 262 4276.6; 263 4267.5; 264 4259.2; 265 4251.9; ...
      266 4245.2; 267 4239.3; 268 4233.9; 269 4229.1; 270 4224.8; ...
      271 4220.9; 272 4217.4; 273 4214.2; 274 4211.3; 275 4208.7; ...
      276 4206.3; 277 4204.1; 278 4202.1; 279 4200.3; 280 4198.7; ...
      281 4197.2; 282 4195.8; 283 4194.5; 284 4193.3; 285 4192.1; ...
      286 4191.1; 287 4190.1; 288 4189.2; 289 4188.4; 290 4187.6; ...
      291 4186.8; 292 4186.1; 293 4185.4; 294 4184.8; 295 4184.2; ...
      296 4183.6; 297 4183.1; 298 4182.6; 299 4182.2; 300 4181.8; ...
      301 4181.4; 302   4181; 303 4180.7; 304 4180.3; 305 4180.1; ...
      306 4179.8; 307 4179.6; 308 4179.4; 309 4179.3; 310 4179.1; ...
      311 4179.1; 312   4179; 313   4179; 314 4178.9; 315   4179; ...
      316   4179; 317 4179.1; 318 4179.2; 319 4179.4; 320 4179.6; ...
      321 4179.8; 322   4180; 323 4180.3; 324 4180.6; 325 4180.9; ...
      326 4181.2; 327 4181.6; 328   4182; 329 4182.5; 330 4182.9; ...
      331 4183.4; 332 4183.9; 333 4184.4; 334 4184.9; 335 4185.5; ...
      336 4186.1; 337 4186.7; 338 4187.3; 339 4187.9; 340 4188.6; ...
      341 4189.2; 342 4189.9; 343 4190.6; 344 4191.2; 345 4191.9; ...
      346 4192.6; 347 4193.3; 348   4194; 349 4194.7; 350 4195.4; ...
      351 4196.1; 352 4196.8; 353 4197.4; 354 4198.1; 355 4198.8; ...
      356 4199.4; 357   4200; 358 4200.6; 359 4201.2; 360 4201.8];

 Tk_water = A(:,1);
 Cp_water = A(:,2);

% interpolate to find Cp_water

 cp = interp1(Tk_water,Cp_water,Tk);

% deal with temperatures beyond range of validity

 if any(Tk < 236)
   L = Tk < 236;

   switch dopt
   case 1
     disp(' ')
     disp('C_water: temperature is beyond range of validity.')
     Tk(L)
     pause
   case 2
     cp(L) = 5699;
   end
 end

 if any(Tk > 360);
   L     = Tk > 360;
   cp(L) = interp1(Tk_water,Cp_water,360);
 end

% find the volumetric heat capacity

 rhocp = rhow * cp;
