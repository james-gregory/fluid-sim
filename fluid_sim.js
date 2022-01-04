PIXI.utils.sayHello();

// N x N grid of cells
let N = 32;
// we will store everything in 1-D arrays
let size = (N + 2) * (N + 2);
// window size
let windowSize = 256;

// grid size
let h = windowSize/(N+2)

// initial arrays
let u = new Array(size).fill(0);          // velocity in x-direction
let v = new Array(size).fill(0);          // velocity in y-direction
let uPrev = new Array(size).fill(0);      // previous u
let vPrev = new Array(size).fill(0);      // previous v
let dens = new Array(size).fill(0);       // density
let densPrev = new Array(size).fill(0);   // previous density
let accel = new Array(size).fill(0);      // acceleration
let accelPrev = new Array(size).fill(0);  // previous acceleration
let rectArray = new Array(size);          // array for storing rectangles


// initial parameters
let strength = 0.1;   // strength of the mouse
let source = 10.0;    // amount of dye added
let visc = 0;         // viscosity
let diff = 0.001;     // diffusivity
let base_u = -0.1;    // base flow

// initalise mouse functions
let mousePositionX = 0.0;
let mousePositionY = 0.0;
let mousePrevX = 0.0;
let mousePrevY = 0.0;
let mouseDown = false;
let rightDown = false;

// function to get the index of (i,j)-th cell
function IX(i_, j_) {
  return i_ + (N+2)*j_;
}

// function to add a source
function addSource(N_, x_, s_, dt_){
  var size = (N_+2) * (N_+2);
  for (let i = 0; i < size; i++){
    x_[i] += dt_*s_[i];
  }
}

// diffuse step
function diffuse(N_, b_, x_, x0_, diff_, dt_){
  a = dt_ * diff_ * N_ * N_;
  for (let k = 0; k < 20; k++){
    for (let i = 1; i <= N_; i++){
      for (let j = 1; j <= N_; j++){
        x_[IX(i, j)] = (x0_[IX(i, j)] + a*(x_[IX(i-1, j)] + x_[IX(i+1, j)] +
                                         x_[IX(i, j-1)] + x_[IX(i, j+1)]))/(1 + 4*a);
      }
    }
    set_bnd(N_, b_, x_);
  }

}

// advect step
function advect(N_, b_, d_, d0_, u_, v_, dt_){
  let dt0 = dt_*N_;
  for (let i = 1; i <= N_; i++){
    for (let j = 1; j <= N_; j++){
      let x = i - dt0*u_[IX(i, j)];
      let y = j - dt0*v_[IX(i, j)];
      if (x < 0.5){x = 0.5;}
      if (x > N_ + 0.5){x = N_ + 0.5;}
      let i0 = Math.floor(x);
      let i1 = i0 + 1;
      if (y < 0.5){y = 0.5;}
      if (y > N_ + 0.5){y = N_ + 0.5;}
      let j0 = Math.floor(y);
      let j1 = j0 + 1;
      let s1 = x - i0;
      let s0 = 1 - s1;
      let t1 = y - j0;
      let t0 = 1 - t1;

      d_[IX(i,j)] = s0*(t0*d0_[IX(i0, j0)] + t1*d0_[i0, j1]) +
                   s1*(t0*d0_[IX(i1, j0)] + t1*d0_[i1, j1]);

    }
  }
  set_bnd(N_, b_, d_);
}

// update velocity step
function updateVelocity(N_, accel_, u_, v_, dt_){

  for (let i = 1; i <= N_; i++){
    for (let j = 1; j <= N_; j++){
      v_[IX(i,j)] += accel[IX(i,j)]*dt_;
    }
  }

  set_bnd(N_, 1, u_);
  set_bnd(N_, 2, v_);

}


// perform a density step
function densStep(N_, x_, x0_, u_, v_, diff_, dt_){
  addSource(N_, x_, x0_, dt_);
  let tmp;

  tmp = x0_;
  x0_ = x_;
  x_ = tmp;

  diffuse(N_, 0, x_, x0_, diff_, dt_);

  tmp = x0_;
  x0_ = x_;
  x_ = tmp;

  advect(N_, 0, x_, x0_, u_, v_, dt_);

}

// perform a velocity step
function velStep(N_, u_, v_, u0_, v0_, visc_, dt_){
  addSource(N_, u_, u0_, dt_);
  addSource(N_, v_, v0_, dt_);

  let tmp;

  tmp = u0_;
  u0_ = u_;
  u_ = tmp;

  diffuse(N_, 1, u_, u0_, visc_, dt_);

  tmp = v0_;
  v0_ = v_;
  v_ = tmp;

  diffuse(N_, 2, v_, v0_, visc_, dt_);
  project(N_, u_, v_, u0_, v0_);


  tmp = u0_;
  u0_ = u_;
  u_ = tmp;

  tmp = v0_;
  v0_ = v_;
  v_ = tmp;

  advect(N_, 2, v_, v0_, u0_, v0_, dt_);
  advect(N_, 1, u_, u0_, u0_, v0_, dt_);
  project(N_, u_, v_, u0_, v0_);
}


function project(N_, u_, v_, p_, div_){
  let h_ = windowSize/(N+2);

  for (let i = 1; i <= N_; i++){
    for (let j = 1; j <= N_; j++){
      div_[IX(i, j)] = -0.5*h_*(u_[IX(i+1, j)] - u_[IX(i-1, j)] + v_[IX(i, j+1)] - v_[IX(i, j-1)]);
      p_[IX(i, j)] = 0;
    }
  }
  set_bnd(N_, 0, div_);
  set_bnd(N_, 0, p_);

  for (let k = 0; k < 20; k++){
    for (let i = 1; i <= N_; i++){
      for (j = 1; j <= N_; j++){
        p_[IX(i, j)] = (div_[IX(i, j)] + p_[IX(i-1, j)] + p_[IX(i+1, j)] +
                                       p_[IX(i, j-1)] + p_[IX(i, j+1)])/4.0;

      }
    }
    set_bnd(N_, 0, p_);
  }

  for (let i = 0; i <= N_; i++){
    for (let j = 1; j <= N_; j++){
      u_[IX(i, j)] -= 0.5*(p_[IX(i+1, j)] - p_[IX(i-1, j)])/h_;
      v_[IX(i, j)] -= 0.5*(p_[IX(i, j+1)] - p_[IX(i, j-1)])/h_;
    }
  }
  set_bnd(N_, 1, u_);
  set_bnd(N_, 2, v_);
}


// set the boundary conditions
// we have the x-velocity equal to base_u at the boundaries and change the
// value of base_u from zero to non-zero
function set_bnd(N_, b_, x_){
  for (let i = 1; i <= N_; i++){
    if (b_ == 0)
    {
      x_[IX(0, i)] = 0;
      x_[IX(N_+1, i)] =0;
      x_[IX(i, 0)] = 0;
      x_[IX(i, N_+1)] = 0;
    }
    else if (b_ == 1)
    {
      x_[IX(0, i)] = base_u;
      x_[IX(N_+1, i)] = base_u;
      x_[IX(i, 0)] = base_u;
      x_[IX(i, N_+1)] = base_u ;
    }
    else if (b_ == 2)
    {
      x_[IX(0, i)] = 0;
      x_[IX(N_+1, i)] =0;
      x_[IX(i, 0)] = 0;
      x_[IX(i, N_+1)] = 0;
    }
  }

  x_[IX(0, 0)] = 0.5 * (x_[IX(1, 0)] + x_[IX(0, 1)]);
  x_[IX(0, N_+1)] = 0.5 * (x_[IX(1, N_+1)] + x_[IX(0, N_)]);
  x_[IX(N_+1, 0)] = 0.5 * (x_[IX(N_, 0)] + x_[IX(N_+1, 1)]);
  x_[IX(N_+1, N_+1)] = 0.5 * (x_[IX(N_, N_+1)] + x_[IX(N_+1, N_)]);
}

// function to change the density based on mouse input
function getFromUI(densPrev_, uPrev_, vPrev_, mouseDown_, rightDown_){

  densPrev_.fill(0);
  uPrev_.fill(0);
  vPrev_.fill(0);

  // do nothing
  if (!mouseDown_ && !rightDown_){
    return;
  }

  // find the cell the mouse is in
  let i = Math.floor((mousePositionX / windowSize) * N + 1);
  let j = Math.floor(((windowSize - mousePositionY) / windowSize) * N + 1);

  // if outside the array do nothing
  if (i < 1 || i > N || j < 1 || j > N){
    return;
  }

  // if left mouse down, add a source AND change the velocity field
  if (mouseDown_){
    densPrev_[IX(i, j)] = source;
    uPrev_[IX(i, j)] = strength * (mousePositionX - mousePrevX);
    vPrev_[IX(i, j)] = strength * -1.0 * (mousePrevY - mousePositionY);
  }

  // if right mouse down, just add a source
  if (rightDown_){
    densPrev_[IX(i, j)] = source;
  }

  mousePrevX = mousePositionX;
  mousePrevY = mousePositionY;

}

function onClick (event) {
  mouseDown = true;
}

function onRelease (event) {
  mouseDown = false;
}

function onMove (event) {
  mousePositionX = event.data.global.x;
  mousePositionY = event.data.global.y;
}

function onRightClick (event) {
  rightDown = true;
}

function onRightRelease (event) {
  rightDown = false;
}

// draw the array on the screen
function drawDens(N_, dens_, densPrev_, tol_, container_){
    for (let i = 1; i < N_+1; i++){
      let x = i * h;
      for (let j = 1; j < N_ + 1; j++){
        let y = (N_ + 1 - j) * h;

          // use the density as the alpha value
          // restrict to between 0 and 1
          var a;
          if (dens_[IX(i, j)] > 1.0)
          {
            a = 1.0;
          }
          else if (dens_[IX(i, j)] < 0.0)
          {
            a = 0.0;
          }
          else
          {
            // the function mapping the density to the alpha value doesn't have
            // to just be alpha = density. We can use different functions.
            //a = dens_[IX(i, j)]**(1);
            // average of surrounding four cells (multiplied by 4!)
            // no particular reason, just settled on through experimentation
            a = (dens_[IX(i, j)] + dens_[IX(i+1, j)]  + dens_[IX(i, j+1)]  + dens_[IX(i+1, j+1)]);
          }

          // store rectangles in an array. Wipe each time we want to draw a new frame
          if (rectArray[IX(i,j)]!= null) rectArray[IX(i,j)].destroy();
          let cell = new PIXI.Graphics();
          rectArray[IX(i,j)] = cell;

          // add a blur filter so we don't just see square cells
          rectArray[IX(i,j)].filters = [blurFilter];

          container_.addChild(rectArray[IX(i,j)]);
          rectArray[IX(i,j)].beginFill(0x000000, a);
          rectArray[IX(i,j)].drawRect(x, y, h, h);
          rectArray[IX(i,j)].endFill();
          container_.addChild(rectArray[IX(i,j)]);

          // debug code for printing the direction vector

          // if (u[IX(i,j)]**2 + u[IX(i,j)]**2 > 0){
          //
          // //console.log(theta);
          // cell.filters = null;
          // [startX, startY, endX, endY] = drawLine(x + 0.5*h, y + 0.5*h, u[IX(i,j)], v[IX(i,j)])
          // cell.lineStyle(2, 0xFF0000);
          // // rectangle.moveTo(x - h*Math.sin(theta), y + h*Math.cos(theta));
          // // rectangle.lineTo(x + h*Math.cos(theta), y + h*Math.sin(theta));
          // cell.moveTo(startX, startY);
          // cell.lineTo(endX, endY);
          //
          // cell.closePath();
          // container_.addChild(cell);
          // }

        //}
      }
    }
}

// function used when drawing the direction vectors on top of the cells
function drawLine(x_, y_, u_, v_){

  var length = h;
  var angle = Math.atan(v_/u_);

  var startX = x_ - 0.5*length*Math.cos(angle);
  var startY = y_ + 0.5*length*Math.sin(angle);

  var endX = x_ + 0.5*length*Math.cos(angle);
  var endY = y_ - 0.5*length*Math.sin(angle);


  return [startX, startY, endX, endY]

}

// function for changing the base flow
function onButtonDown()
{
  if (base_u != 0){
    base_u = 0;
  }
  else base_u = -0.1;
}



// Create the application helper and add its render target to the page
let app = new PIXI.Application({ width: windowSize+96, height: windowSize , backgroundColor: 0xFFFFFF});
document.body.appendChild(app.view);

let elapsed = 0.0;

// load in the fan image and turn into a button
var textureButton = PIXI.Texture.from('fan.jpg');
var button = new PIXI.Sprite(textureButton);
   button.buttonMode = true;

   button.anchor.set(0.5);

   button.position.x = windowSize + 25;
   button.position.y = windowSize*0.5;

   // make the button interactive...
   button.interactive = true;
   button.on('mousedown', onButtonDown)

app.stage.addChild(button);

// define the blur filter (gaussian)
var blurFilter = new PIXI.filters.BlurFilter();
blurFilter.blur = 10;
blurFilter.quality = 7;


//ticker to update app
app.ticker.maxFPS = 60;
app.ticker.speed = 0.1;
app.ticker.add((delta) => {

  // update time
  elapsed += delta;

  let container = new PIXI.Container();
  app.renderer.plugins.interaction.on('mousemove', onMove);
  app.renderer.plugins.interaction.on('mousedown', onClick);
  app.renderer.plugins.interaction.on('mouseup', onRelease);
  app.renderer.plugins.interaction.on('rightdown', onRightClick);
  app.renderer.plugins.interaction.on('rightup', onRightRelease);
  getFromUI(densPrev, uPrev, vPrev, mouseDown, rightDown);

  velStep(N, u, v, uPrev, vPrev, visc, delta);
  densStep(N, dens, densPrev, u, v, diff, delta);
  drawDens(N, dens, densPrev, 0.01, container);
  app.stage.addChild(container);

});
