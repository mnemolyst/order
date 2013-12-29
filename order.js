var drawText = false;
var renderInterval = 33; //milliseconds
var gravity = false;
var gravConst = 32.2 * 12 / 1000000 / 8.25; // ft/s^2 -> heads/ms^2 (1 head = 8.25 in)
var framerate = 0;
var hoist = false;
var pause = false;

function crossProd(x, y) {
  return [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]];
}

function Scene2D(canvas) {
  this.matrix = new Matrix2D();
  this.matrixGood = false;
  this.context = canvas.getContext('2d');
  this.sceneWidth = canvas.width;
  this.sceneHeight = canvas.height;
  this.panX = 0;
  this.panY = 0;
  this.scale = 1.0;
  this.rotation = 0;
  this.points = [];
  this.pointsTransformed = [];
  this.pointsLast = [];
  this.pointsAccel = [];
  this.pointsInvMass = [];
  this.numPoints = 0;
  this.polyItems = [];
  this.numPolyItems = 0;
  this.constraints = [];
  this.numConstraints = 0;
  this.speed = 0.1;
  this.constrainOrder = [];
  this.dragItem = null;
  this.deleteDrag = false;

  this.addPoint = function(x, y) {
    this.points[this.points.length] = this.pointsLast[this.pointsLast.length] = x;
    this.points[this.points.length] = this.pointsLast[this.pointsLast.length] = y;
    this.pointsAccel[this.pointsAccel.length] = 0;
    this.pointsAccel[this.pointsAccel.length] = 0;
    this.pointsInvMass[this.pointsInvMass.length] = 1;

    return this.numPoints++;
  }

  this.getPoint = function(idx) {
    var idx2D = idx*2;
    return [this.points[idx2D], this.points[idx2D+1]];
  }

  this.addPolyItem = function(item) {
    this.polyItems[this.polyItems.length] = item;

    return this.numPolyItems++;
  }

  this.randomizeConstrainOrder = function() {
    for (var i=0; i<this.polyItems.length; i++) {
      this.constrainOrder[i] = i;
    }

    for (var i=this.constrainOrder.length-1; i>0; i--) {
      var j = Math.floor(Math.random() * (i+1));
      var temp = this.constrainOrder[i];
      this.constrainOrder[i] = this.constrainOrder[j];
      this.constrainOrder[j] = temp;
    }
  }

  this.render = function() {
    var halfWidth = this.sceneWidth*0.5;
    var halfHeight = this.sceneHeight*0.5;

    if (!this.matrixGood) {
      this.matrix.identity();
      this.matrix.translate(halfWidth-this.panX, halfHeight-this.panY);
      this.matrix.scale(this.scale, this.scale);
      //this.matrix.rotate(this.rotation);
      this.matrixGood = true;
    }

    this.pointsTransformed = this.matrix.transformArray(this.points);

    this.context.save();
    this.context.fillStyle = '#000000';
    this.context.fillRect(0, 0, this.sceneWidth, this.sceneHeight);
    if (drawText) {
      this.context.fillStyle = '#aaaaaa';
      this.context.textBaseline = 'top';
      this.context.fillText('framerate: ' + framerate, 5, 5);
      this.context.fillText('polys: ' + this.numPolyItems, 5, 20);
    }

    for (var i=0; i<this.numPolyItems; i++) {
      this.polyItems[i].render(this);
    }

    //if (this.dragItem != null)
    //  this.dragItem.render(this);

    this.context.restore();
  }

  this.accumulateForces = function() {
    for (var i=0; i<this.numPoints; i++) {
      var i2 = i*2;
      var xIdx = i2;
      var yIdx = i2+1;

      //this.points3DAccel[xIdx] = gravity?gravConst*windCos*0.1:0;
      //this.points3DAccel[yIdx] = gravity?gravConst*windSin*0.1:0;
      this.pointsAccel[xIdx] = 0;
      this.pointsAccel[yIdx] = gravity?-gravConst:0; // 1g of gravity in heads (8.25in: female) per millisec^2
    }

    this.snapPolys();
  }

  this.verlet = function(timeDel) {
    for (var i=0; i<this.numPoints; i++) {
      var i2 = i*2;
      var xIdx = i2;
      var yIdx = i2+1;

      var tempX = this.points[xIdx];
      var tempY = this.points[yIdx];

      this.points[xIdx] += (this.points[xIdx]-this.pointsLast[xIdx])*0.99 + this.pointsAccel[xIdx]*timeDel*timeDel;
      this.points[yIdx] += (this.points[yIdx]-this.pointsLast[yIdx])*0.99 + this.pointsAccel[yIdx]*timeDel*timeDel;
      //this.points[xIdx] += this.points[xIdx]-this.pointsLast[xIdx] + this.pointsAccel[xIdx]*timeDel*timeDel;
      //this.points[yIdx] += this.points[yIdx]-this.pointsLast[yIdx] + this.pointsAccel[yIdx]*timeDel*timeDel;

      this.pointsLast[xIdx] = tempX;
      this.pointsLast[yIdx] = tempY;
    }
  }

  this.constrain = function() {
    for (var i=0; i<this.numPolyItems; i++) {
      this.polyItems[i].constrain();
    }
    if (this.dragItem != null)
      this.dragItem.constrain();
  }

  this.clip = function() {
    var temp = 0;
    for (var i=0; i<this.numPoints; i++) {
      if (this.points[i*2] < -this.sceneWidth*0.5 || this.points[i*2] > this.sceneWidth*0.5) {
        temp = this.points[i*2];
        this.points[i*2] = this.pointsLast[i*2];
        this.pointsLast[i*2] = temp;
      }
      if (this.points[i*2+1] < -this.sceneHeight*0.5 || this.points[i*2+1] > this.sceneHeight*0.5) {
        temp = this.points[i*2+1];
        this.points[i*2+1] = this.pointsLast[i*2+1];
        this.pointsLast[i*2+1] = temp;
      }
    }
  }

  this.collidePolys = function() {
    //for (var i=0; i<this.numPolyItems; i++) {
    //  this.polyItems[i].color = '#eeeeee';
    //}
    for (var i=0; i<this.numPolyItems; i++) {
      for (var j=i+1; j<this.numPolyItems; j++) {
        this.polyItems[i].collidePoly(this.polyItems[j]);
        //if (this.polyItems[i].collidePoly(this.polyItems[j])) {
        //  this.polyItems[i].color = '#ee0000';
        //  this.polyItems[j].color = '#ee0000';
        //}
      }
    }
  }

  this.snapPolys = function() {
    for (var i=0; i<this.numPolyItems; i++) {
      this.polyItems[i].color = '#eeeeee';
    }
    for (var i=0; i<this.numPolyItems; i++) {
      for (var j=i+1; j<this.numPolyItems; j++) {
        //this.polyItems[i].snapPoly(this.polyItems[j]);
        if (this.polyItems[i].snapPoly(this.polyItems[j])) {
          this.polyItems[i].color = '#ee0000';
          this.polyItems[j].color = '#ee0000';
        }
      }
    }
  }

  this.safeActions = function() {
    if (this.deleteDrag) {
      this.dragItem = null;
      this.deleteDrag = false;
    }
  }

  this.tic = function() {
    this.render();
    this.accumulateForces();
    this.verlet(renderInterval * this.speed);
    this.constrain();
    this.collidePolys();
    this.collidePolys();
    this.constrain();
    this.collidePolys();
    this.collidePolys();
    this.clip();
    this.safeActions();
  }
}

function polyItem(scene, pIdx) {
  this.scene = scene;
  this.pIdx = pIdx;
  this.p2Didx = [];
  for (var i=0; i<this.pIdx.length; i++) {
    this.p2Didx[i] = 2 * this.pIdx[i];
  }
  this.sideLengths = [];
  this.invSideLengths = [];
  this.restLengths = [];
  this.restLengthsSq = [];
  this.edgeVecIdx = [];
  this.normalDepths = [];
  this.color = '#eeeeee';

  // INIT
  var points = this.scene.points;
  var minSelfProj = 0;
  for (var i=0; i<this.p2Didx.length; i++) {
    var si = (i+1<this.p2Didx.length)?i+1:0; // successor of i

    // compute side lengths
    var v = [points[this.p2Didx[si]]-points[this.p2Didx[i]], points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1]];
    this.sideLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
    this.invSideLengths[i] = 1.0 / this.sideLengths[i];

    // constraint length pairs
    this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[i];
    this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[si];
    var n = this.p2Didx.length-3-(i>1?i-1:0);
    for (var j=i+2; j<i+n+2; j++) {
      this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[i];
      this.edgeVecIdx[this.edgeVecIdx.length] = this.p2Didx[j];
    }

    // depth from each face
    var norm = [(points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i], // right-hand unit normal of i-th side
                (points[this.p2Didx[i]]-points[this.p2Didx[si]])*this.invSideLengths[i]];
    this.normalDepths[i] = 0;
    for (var j=0; j<this.p2Didx.length; j++) {
      var proj = (points[this.p2Didx[j]]-points[this.p2Didx[i]])*norm[0] + (points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])*norm[1];
      this.normalDepths[i] = (this.normalDepths[i] < -proj)?-proj:this.normalDepths[i];
    }
  }

  // constraint rest lengths
  for (var i=0; i<this.edgeVecIdx.length/2; i++) {
    var v = [points[this.edgeVecIdx[i*2+1]]-points[this.edgeVecIdx[i*2]], points[this.edgeVecIdx[i*2+1]+1]-points[this.edgeVecIdx[i*2]+1]];
    this.restLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
    this.restLengthsSq[i] = this.restLengths[i]*this.restLengths[i];
  }

  // METHODS
  this.render = function() {
    var points = this.scene.pointsTransformed;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = this.color;
    for (var i=0; i<this.p2Didx.length; i++) {
      context.lineTo(points[this.p2Didx[i]], this.scene.sceneHeight - points[this.p2Didx[i]+1]);
    }
    context.lineTo(points[this.p2Didx[0]], this.scene.sceneHeight - points[this.p2Didx[0]+1]); // close the polygon
    context.stroke();
  }

  this.constrain = function() {
    var points = this.scene.points;
    for (var i=0; i<this.edgeVecIdx.length/2; i++) { //TODO 
      var vecH = i*2;
      var delta = [points[this.edgeVecIdx[vecH+1]]-points[this.edgeVecIdx[vecH]],
                   points[this.edgeVecIdx[vecH+1]+1]-points[this.edgeVecIdx[vecH]+1]];
      var delSq = delta[0]*delta[0] + delta[1]*delta[1];
      var im1 = 1;//this.scene.pointsInvMass[this.edgeVecIdx[vecH]/2];
      var im2 = 1;//this.scene.pointsInvMass[this.edgeVecIdx[vecH+1]/2];
      var diff = (delSq-this.restLengthsSq[i])/((delSq+this.restLengthsSq[i])*(im1+im2)*2);
      for (var j=0; j<2; j++) {
        var del = diff*delta[j];
        points[this.edgeVecIdx[vecH]+j]   += del*im1;
        points[this.edgeVecIdx[vecH+1]+j] -= del*im2;
      }
    }
  }

  this.collidePoly = function(poly2) {
    var points = this.scene.points;

    var minOverlap = 1000;
    var bestNorm = [];
    var besti = -1;
    var bestj = -1;
    var bestp = -1;

    // First check this poly's sides for penetration against "that" poly's points
    for (var i=0; i<this.p2Didx.length; i++) {
      var j = (i+1<this.p2Didx.length)?i+1:0;
      var norm = [(points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i], // unit normal of i-th side
                  (points[this.p2Didx[i]]-points[this.p2Didx[j]])*this.invSideLengths[i]];

      var min12Proj = 0;
      var max12Proj = 0;
      var min12Projp = -1;
      for (var p=0; p<poly2.p2Didx.length; p++) {
        var proj = (points[poly2.p2Didx[p]]-points[this.p2Didx[i]])*norm[0] + (points[poly2.p2Didx[p]+1]-points[this.p2Didx[i]+1])*norm[1];
        if (proj<min12Proj) {
          min12Proj = proj;
          min12Projp = p;
        }
        max12Proj = (max12Proj<proj)?proj:max12Proj;
      }

      var minProj = (-this.normalDepths[i]<min12Proj)?-this.normalDepths[i]:min12Proj;
      var maxProj = (max12Proj>0)?max12Proj:0;

      if (min12Projp == -1)
        return false;

      var overlap = max12Proj - min12Proj + this.normalDepths[i] - maxProj + minProj;
      if (overlap > 0) {
        //alert(overlap);
        if (overlap < minOverlap) {
          minOverlap = overlap;
          bestNorm = norm;
          besti = this.p2Didx[i];
          bestj = this.p2Didx[j];
          bestp = poly2.p2Didx[min12Projp];
        }
      } else return false;
    }

    // Then check "that" poly's sides
    for (var i=0; i<poly2.p2Didx.length; i++) {
      var j = (i+1<poly2.p2Didx.length)?i+1:0;
      var norm = [(points[poly2.p2Didx[j]+1]-points[poly2.p2Didx[i]+1])*poly2.invSideLengths[i], // unit normal of i-th side
                  (points[poly2.p2Didx[i]]-points[poly2.p2Didx[j]])*poly2.invSideLengths[i]];

      var min21Proj = 0;
      var max21Proj = 0;
      var min21Projp = -1;
      for (var p=0; p<this.p2Didx.length; p++) {
        var proj = (points[this.p2Didx[p]]-points[poly2.p2Didx[i]])*norm[0] + (points[this.p2Didx[p]+1]-points[poly2.p2Didx[i]+1])*norm[1];
        if (proj<min21Proj) {
          min21Proj = proj;
          min21Projp = p;
        }
        max21Proj = (max21Proj<proj)?proj:max21Proj;
      }

      var minProj = (-poly2.normalDepths[i]<min21Proj)?-poly2.normalDepths[i]:min21Proj;
      var maxProj = (max21Proj>0)?max21Proj:0;

      if (min21Projp == -1)
        return false;

      var overlap = max21Proj - min21Proj + poly2.normalDepths[i] - maxProj + minProj;
      if (overlap > 0) {
        //alert(overlap);
        if (overlap < minOverlap) {
          minOverlap = overlap;
          bestNorm = norm;
          besti = poly2.p2Didx[i];
          bestj = poly2.p2Didx[j];
          bestp = this.p2Didx[min21Projp];
        }
      } else return false;
    }

    // Update points
    if (Math.abs(points[bestj] - points[besti]) > Math.abs(points[bestj+1] - points[besti+1]))
      var ijHitF = Math.abs(points[bestp]-points[besti]) / Math.abs(points[bestj]-points[besti]);
    else
      var ijHitF = Math.abs(points[bestp+1]-points[besti+1]) / Math.abs(points[bestj+1]-points[besti+1]);
    var l = 1.0/(ijHitF*ijHitF + (1-ijHitF)*(1-ijHitF));
    var iHitP = (1-ijHitF)*l;
    var jHitP = ijHitF*l;
    var push = [bestNorm[0]*minOverlap*0.5, bestNorm[1]*minOverlap*0.5];
    //alert (bestNorm[0] + ' ' + bestNorm[1]);
    points[besti]   -= push[0]*iHitP;
    points[besti+1] -= push[1]*iHitP;
    points[bestj]   -= push[0]*jHitP;
    points[bestj+1] -= push[1]*jHitP;
    points[bestp]   += push[0];
    points[bestp+1] += push[1];

    return true;
  }

  this.snapPoly = function(poly2) {
    var points = this.scene.points;
    var pointsLast = this.scene.pointsLast;
    var pointsAccel = this.scene.pointsAccel;

    var minOverlap = 1000;
    var bestNorm = [];
    var besti = -1;
    var bestj = -1;
    var bestp = -1;

    // check all this poly's sides against all that poly's sides
    for (var i=0; i<this.p2Didx.length; i++) {
      var si = (i+1<this.p2Didx.length)?i+1:0;
      for (var j=0; j<poly2.p2Didx.length; j++) {
        var sj = (j+1<poly2.p2Didx.length)?j+1:0;

        var v1 = [points[poly2.p2Didx[sj]]-points[this.p2Didx[i]], points[poly2.p2Didx[sj]+1]-points[this.p2Didx[i]+1]];
        var v2 = [points[poly2.p2Didx[j]]-points[this.p2Didx[si]], points[poly2.p2Didx[j]+1]-points[this.p2Didx[si]+1]];
        var ds1 = v1[0]*v1[0]+v1[1]*v1[1]; // square of distance from this poly's point i to that poly's point sj
        var ds2 = v2[0]*v2[0]+v2[1]*v2[1];
        if (ds1<400.0 && ds2<400.0) {
          var d1 = Math.sqrt(ds1);
          var d2 = Math.sqrt(ds2);
          var n1 = [v1[1], -v1[0]];
          var n2 = [v2[1], -v2[0]];
          var ui = [(points[this.p2Didx[si]]-points[this.p2Didx[i]])*this.invSideLengths[i], (points[this.p2Didx[si]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i]];
          var proj1 = v1[0]*ui[0] + v1[1]*ui[1];
          var proj2 = v2[0]*ui[0] + v2[1]*ui[1];

          pointsAccel[this.p2Didx[i]] +=     (v1[0]*(d1-10.0) - n1[0]*proj1)*0.001 -  (points[this.p2Didx[i]] - pointsLast[this.p2Didx[i]])*0.01;
          pointsAccel[this.p2Didx[i]+1] +=   (v1[1]*(d1-10.0) - n1[1]*proj1)*0.001 -  (points[this.p2Didx[i]+1] - pointsLast[this.p2Didx[i]+1])*0.01;
          pointsAccel[this.p2Didx[si]] +=    (v2[0]*(d2-10.0) - n2[0]*proj2)*0.001 -  (points[this.p2Didx[si]] - pointsLast[this.p2Didx[si]])*0.01;
          pointsAccel[this.p2Didx[si]+1] +=  (v2[1]*(d2-10.0) - n2[1]*proj2)*0.001 -  (points[this.p2Didx[si]+1] - pointsLast[this.p2Didx[si]+1])*0.01;
          pointsAccel[poly2.p2Didx[j]] +=    (-v2[0]*(d2-10.0) + n2[0]*proj2)*0.001 - (points[poly2.p2Didx[j]] - pointsLast[poly2.p2Didx[j]])*0.01;
          pointsAccel[poly2.p2Didx[j]+1] +=  (-v2[1]*(d2-10.0) + n2[1]*proj2)*0.001 - (points[poly2.p2Didx[j]+1] - pointsLast[poly2.p2Didx[j]+1])*0.01;
          pointsAccel[poly2.p2Didx[sj]] +=   (-v1[0]*(d1-10.0) + n1[0]*proj1)*0.001 - (points[poly2.p2Didx[sj]] - pointsLast[poly2.p2Didx[sj]])*0.01;
          pointsAccel[poly2.p2Didx[sj]+1] += (-v1[1]*(d1-10.0) + n1[1]*proj1)*0.001 - (points[poly2.p2Didx[sj]+1] - pointsLast[poly2.p2Didx[sj]+1])*0.01;

          return true;
        }
      }
    }

    return false;
  }

  this.pointIntersects = function(point) {
    var points = this.scene.points;

    for (var i=0; i<this.p2Didx.length; i++) {
      var j = (i+1<this.p2Didx.length)?i+1:0;
      var norm = [(points[this.p2Didx[j]+1]-points[this.p2Didx[i]+1])*this.invSideLengths[i], // unit normal of i-th side
                  (points[this.p2Didx[i]]-points[this.p2Didx[j]])*this.invSideLengths[i]];
      var proj = (point[0]-points[this.p2Didx[i]])*norm[0] + (point[1]-points[this.p2Didx[i]+1])*norm[1];
      if (proj > 0)
        return false;
    }

    return true;
  }
}

function dragItem(scene, dragP, pIdx) {
  this.scene = scene;
  this.dragP = dragP;
  this.pIdx = pIdx;
  this.p2Didx = [];
  for (var i=0; i<this.pIdx.length; i++) {
    this.p2Didx[i] = 2 * this.pIdx[i];
  }
  this.restLengths = [];
  this.restLengthsSq = [];
  this.color = '#eeeeee';

  // INIT
  var points = this.scene.points;
  // constraint rest lengths
  for (var i=0; i<this.p2Didx.length; i++) {
    var v = [points[this.p2Didx[i]]-this.dragP[0], points[this.p2Didx[i]+1]-this.dragP[1]];
    this.restLengths[i] = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
    this.restLengthsSq[i] = this.restLengths[i]*this.restLengths[i];
  }

  // METHODS
  this.render = function() {
    var points = this.scene.pointsTransformed;
    var context = this.scene.context;
    var dragPointTransformed = this.scene.matrix.transformArray(this.dragP);

    context.beginPath();
    context.strokeStyle = this.color;
    for (var i=0; i<this.p2Didx.length; i++) {
      context.moveTo(dragPointTransformed[0], this.scene.sceneHeight - dragPointTransformed[1]);
      context.lineTo(points[this.p2Didx[i]], this.scene.sceneHeight - points[this.p2Didx[i]+1]);
    }
    context.stroke();
  }

  this.constrain = function() {
    var points = this.scene.points;
    for (var i=0; i<this.p2Didx.length; i++) {
      var delta = [points[this.p2Didx[i]]-this.dragP[0],
                   points[this.p2Didx[i]+1]-this.dragP[1]]
      var delSq = delta[0]*delta[0] + delta[1]*delta[1];
      var diff = (delSq-this.restLengthsSq[i])/(delSq+this.restLengthsSq[i]);
      for (var j=0; j<2; j++) {
        var del = diff*delta[j];
        points[this.p2Didx[i]+j] -= del;
      }
    }
  }
}

function linItem(scene, p1idx, p2idx) {
  this.scene = scene;
  this.p1idx = p1idx;
  this.p2idx = p2idx;
  this.p12Didx = p1idx*2;
  this.p22Didx = p2idx*2;
  this.p13Didx = p1idx*3;
  this.p23Didx = p2idx*3;

  // INIT
  // compute rest length
  var p1 = scene.getPoint3D(p1idx);
  var p2 = scene.getPoint3D(p2idx);
  this.restLength = Math.sqrt((p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]));
  this.restSq = this.restLength * this.restLength;

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#eeeeee';
    context.lineTo(points[this.p12Didx], points[this.p12Didx+1]);
    context.lineTo(points[this.p22Didx], points[this.p22Didx+1]);
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var delta = [points[this.p23Didx]-points[this.p13Didx], points[this.p23Didx+1]-points[this.p13Didx+1], points[this.p23Didx+2]-points[this.p13Didx+2]];
    var delSq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
    var diff = (delSq-this.restSq)/(delSq+this.restSq)*0.5;
    for (var j=0; j<3; j++) {
      points[this.p13Didx+j] += diff * delta[j];
      points[this.p23Didx+j] -= diff * delta[j];
    }
  }
}

function cordItem(scene, maxLength, points) {
  this.scene = scene;
  this.maxLength = maxLength;
  this.avgMaxSeg = maxLength / (points.length - 1);
  this.points = points;
  this.p2Didx = [];
  this.p3Didx = [];
  for (var i=0; i<points.length; i++) {
    this.p2Didx[i] = points[i]*2;
    this.p3Didx[i] = points[i]*3;
  }

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#ee6666';
    for (var i=0; i<this.p2Didx.length; i++) {
      context.lineTo(points[this.p2Didx[i]], points[this.p2Didx[i]+1]);
    }
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var deltas = [];
    var delSqs = [];
    var delLens = [];
    var totalLen = 0;
    for (var i=1; i<this.p3Didx.length; i++) {
      var delta = deltas[i-1] = [points[this.p3Didx[i]]-points[this.p3Didx[i-1]],
                                 points[this.p3Didx[i]+1]-points[this.p3Didx[i-1]+1],
                                 points[this.p3Didx[i]+2]-points[this.p3Didx[i-1]+2]];
      delSqs[i-1] = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      //totalLen += delLens[i-1] = Math.sqrt(delSqs[i-1]);
      var apx = delSqs[i-1];
      apx = (apx + delSqs[i-1]/apx)*0.5;
      apx = (apx + delSqs[i-1]/apx)*0.5;
      apx = (apx + delSqs[i-1]/apx)*0.5;
      totalLen += delLens[i-1] = apx;
    }
    if (totalLen <= maxLength)
      return;

    var err = (totalLen-maxLength) / (this.p3Didx.length-1);
    for (var i=1; i<this.p3Didx.length; i++) {
      var diff = err/(delLens[i-1]*2);
      for (var j=0; j<3; j++) {
        points[this.p3Didx[i-1]+j] += diff * deltas[i-1][j];
        points[this.p3Didx[i]+j] -= diff * deltas[i-1][j];
      }
    }
  }
}

function linSumPairItem(scene, p1idx, p2idx, p3idx, p4idx, maxShrink, maxGrow) {
  this.scene = scene;
  this.p1idx = p1idx;
  this.p2idx = p2idx;
  this.p3idx = p3idx;
  this.p4idx = p4idx;
  this.maxShrink = maxShrink;
  this.maxGrow = maxGrow;
  this.p12Didx = p1idx*2;
  this.p22Didx = p2idx*2;
  this.p32Didx = p3idx*2;
  this.p42Didx = p4idx*2;
  this.p13Didx = p1idx*3;
  this.p23Didx = p2idx*3;
  this.p33Didx = p3idx*3;
  this.p43Didx = p4idx*3;

  // INIT
  // compute rest length
  var p1 = scene.getPoint3D(p1idx);
  var p2 = scene.getPoint3D(p2idx);
  var p3 = scene.getPoint3D(p3idx);
  var p4 = scene.getPoint3D(p4idx);
  this.restLength = Math.sqrt((p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]));
  this.restLength += Math.sqrt((p4[0]-p3[0])*(p4[0]-p3[0]) + (p4[1]-p3[1])*(p4[1]-p3[1]) + (p4[2]-p3[2])*(p4[2]-p3[2]));
  this.restSq = this.restLength * this.restLength;
  this.minLength = this.restLength - maxShrink;
  this.maxLength = this.restLength + maxGrow;

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#ee3333';
    context.lineTo(points[this.p12Didx], points[this.p12Didx+1]);
    context.lineTo(points[this.p22Didx], points[this.p22Didx+1]);
    context.stroke();
    context.moveTo(points[this.p32Didx], points[this.p32Didx+1]);
    context.lineTo(points[this.p42Didx], points[this.p42Didx+1]);
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var delta1 = [points[this.p23Didx]-points[this.p13Didx], points[this.p23Didx+1]-points[this.p13Didx+1], points[this.p23Didx+2]-points[this.p13Didx+2]];
    var delta2 = [points[this.p43Didx]-points[this.p33Didx], points[this.p43Didx+1]-points[this.p33Didx+1], points[this.p43Didx+2]-points[this.p33Didx+2]];
    var delSq1 = delta1[0]*delta1[0] + delta1[1]*delta1[1] + delta1[2]*delta1[2];
    var delSq2 = delta2[0]*delta2[0] + delta2[1]*delta2[1] + delta2[2]*delta2[2];
    var delLen = Math.sqrt(delSq1) + Math.sqrt(delSq2); //TODO optimize

    if (delLen > this.maxLength) {
      var diff = (delLen - this.maxLength) / delLen * 0.23;
    } else if (delLen < this.minLength) {
      var diff = (delLen - this.minLength) / delLen * 0.23;
    } else {
      return;
    }

    for (var j=0; j<3; j++) {
      points[this.p13Didx+j] += diff * delta1[j];
      points[this.p23Didx+j] -= diff * delta1[j];
      points[this.p33Didx+j] += diff * delta2[j];
      points[this.p43Didx+j] -= diff * delta2[j];
    }
  }
}

function linDiffPairItem(scene, p1idx, p2idx, p3idx, p4idx, maxDiff) {
  this.scene = scene;
  this.p1idx = p1idx;
  this.p2idx = p2idx;
  this.p3idx = p3idx;
  this.p4idx = p4idx;
  this.maxDiff = maxDiff;
  this.p12Didx = p1idx*2;
  this.p22Didx = p2idx*2;
  this.p32Didx = p3idx*2;
  this.p42Didx = p4idx*2;
  this.p13Didx = p1idx*3;
  this.p23Didx = p2idx*3;
  this.p33Didx = p3idx*3;
  this.p43Didx = p4idx*3;

  // INIT
  // compute rest length
  var p1 = scene.getPoint3D(p1idx);
  var p2 = scene.getPoint3D(p2idx);
  var p3 = scene.getPoint3D(p3idx);
  var p4 = scene.getPoint3D(p4idx);
  this.restDiff = Math.sqrt((p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]));
  this.restDiff -= Math.sqrt((p4[0]-p3[0])*(p4[0]-p3[0]) + (p4[1]-p3[1])*(p4[1]-p3[1]) + (p4[2]-p3[2])*(p4[2]-p3[2]));
  this.restDiff = Math.abs(this.restDiff);
  this.restSq = this.restDiff * this.restDiff;

  // METHODS
  this.render = function() {
    var points = this.scene.points2D;
    var context = this.scene.context;

    context.beginPath();
    context.strokeStyle = '#33ee33';
    context.lineTo(points[this.p12Didx], points[this.p12Didx+1]);
    context.lineTo(points[this.p22Didx], points[this.p22Didx+1]);
    context.stroke();
    context.moveTo(points[this.p32Didx], points[this.p32Didx+1]);
    context.lineTo(points[this.p42Didx], points[this.p42Didx+1]);
    context.stroke();
  }

  this.constrain= function() {
    var points = this.scene.points3D;
    var delta1 = [points[this.p23Didx]-points[this.p13Didx], points[this.p23Didx+1]-points[this.p13Didx+1], points[this.p23Didx+2]-points[this.p13Didx+2]];
    var delta2 = [points[this.p43Didx]-points[this.p33Didx], points[this.p43Didx+1]-points[this.p33Didx+1], points[this.p43Didx+2]-points[this.p33Didx+2]];
    var delSq1 = delta1[0]*delta1[0] + delta1[1]*delta1[1] + delta1[2]*delta1[2];
    var delSq2 = delta2[0]*delta2[0] + delta2[1]*delta2[1] + delta2[2]*delta2[2];
    var delLen = Math.sqrt(delSq1) - Math.sqrt(delSq2); //TODO optimize

    if (delLen > this.maxDiff) {
      var diff = (delLen - this.maxDiff) / delLen * 0.23;
    } else if (delLen < -this.maxDiff) {
      var diff = (-delLen - this.maxDiff) / delLen * 0.23;
    } else {
      return;
    }

    for (var j=0; j<3; j++) {
      points[this.p13Didx+j] += diff * delta1[j];
      points[this.p23Didx+j] -= diff * delta1[j];
      points[this.p33Didx+j] -= diff * delta2[j];
      points[this.p43Didx+j] += diff * delta2[j];
    }
  }
}

jQuery(document).ready(function($) {
  var canvas = $('#baseCanvas');
  var offset = canvas.offset();
  var scene = new Scene2D(canvas[0]);
  var state = 'none';
  var px = 0;
  var py = 0;
  var i;

  function addPoly(n) {
    //var x = Math.random() * canvas.width() - canvas.width()/2;
    //var y = Math.random() * canvas.height() - canvas.height()/2;
    var x = px;
    var y = py;
    var r = 50;
    var t = 2*Math.PI/n

    var pIdx = [];
    for (var i=0; i<n; i++) {
      pIdx[i] = scene.addPoint(x+r*Math.cos(t/2 + i*t), y+r*Math.sin(t/2 + i*t));
    }
    scene.addPolyItem(new polyItem(scene, pIdx));
  }

  //var points = {};
  //for (var p in pointData) {
  //  var point = pointData[p];
  //  points[p] = scene.addPoint(point[0], point[1]);
  //}

  //for (var item in polyItemData) {
  //  var pList = polyItemData[item]['points'];
  //  scene.addPolyItem(new polyItem(scene, [points[pList[0]], points[pList[1]], points[pList[2]], points[pList[3]]]));
  //}

  //for (var item in linSumPairItemData) {
  //  var data = linSumPairItemData[item];
  //  scene.addPolyItem(new linSumPairItem(scene, points[data[0]], points[data[1]], points[data[2]], points[data[3]], data[4], data[5]));
  //}

  //for (var item in linDiffPairItemData) {
  //  var data = linDiffPairItemData[item];
  //  scene.addPolyItem(new linDiffPairItem(scene, points[data[0]], points[data[1]], points[data[2]], points[data[3]], data[4]));
  //}

  var frames=0;

  setInterval(function () {
    scene.tic();
    frames++;
  }, renderInterval);

  setInterval(function () {
    if (drawText) {
      framerate = frames;
      frames=0;
    }
  }, 1000);

  //setInterval(function () {
  //  if (scene.numPolyItems < 10)
  //    addQuad();
  //}, 4000);

  $(document).keydown(function(event) {
    if (state == 'alert_next_key') {
      state = 'none';
      alert(event.which);
    }

    // KEYS
    switch(event.which) {
      case 8:                                       //backspace
        break;
      case 9:                                       //tab
        break;
      case 13:                                      //enter
        break;
      case 27:                                      //esc
        break;
      case 37:                                      //left arrow
        if (event.shiftKey) {
        } else {
          scene.panX -= 5;
        }
        break;
      case 38:                                      //up arrow
        scene.panY -= 5;
        break;
      case 39:                                      //right arrow
        if (event.shiftKey) {
        } else {
          scene.panX += 5;
        }
        break;
      case 40:                                      //down arrow
        scene.panY += 5;
        break;
      case 46:                                      //delete
        break;
      case 48: case 49: case 50:                    //'0'-'2'
        break;
      case 51: case 52: case 53: case 54: case 55: case 56: case 57: //'3'-'9'
        addPoly(event.which - 48);
        break;
      case 61:                                      //'=' / '+'
        break;
      case 66:                                      //'b' bomb
        var x = Math.random() * canvas.width() - canvas.width()/2;
        var y = Math.random() * canvas.height() - canvas.height()/2;
        var pow = 0.02;

        for (var i=0; i<scene.points3D.length; i+=3) {
          var delta = [scene.points3D[i]-x, scene.points3D[i+1]-y, scene.points3D[i+2]-z];
          var deltalength = Math.sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]); //TODO optimize
          var c = pow/(deltalength*deltalength*deltalength);
          for (var j=0; j<3; j++) {
            scene.points3D[i+j] += delta[j]*c;
          }
        }
        break;
      case 67:                                      //'c'
        break;
      case 68:                                      //'d'
        addDart(); //TODO
        break;
      case 69:                                      //'e'
        break;
      case 70:                                      //'f'
        addPoly(Math.floor(Math.random()*6+3));
        break;
      case 71:                                      //'g'
        gravity = !gravity;
        break;
      case 72:                                      //'h'
        hoist = !hoist;
        break;
      case 73:                                      //'i'
        scene.scale *= 1.1;
        break;
      case 76:                                      //'l'
        break;
      case 78:                                      //'n'
        drawText = !drawText;
        break;
      case 79:                                      //'o'
        scene.scale /= 1.1;
        break;
      case 80:                                      //'p'
        pause = !pause;
        break;
      case 82:                                      //'r'
        break;
      case 83:                                      //'s'
        break;
      case 85:                                      //'u'
        scene.zTranslate -= 0.2;
        break;
      case 86:                                      //'v'
        break;
      case 87:                                      //'w'
        break;
      case 89:                                      //'y'
        break;
      case 90:                                      //'z'
        break;
      case 173:                                     //'-'
        break;
      case 190:                                     //'.'
        break;
      case 191:                                     //'/'
        break;
      case 192:                                     //'`'
        if (event.shiftKey)
          state = 'alert_next_key';
        break;
      case 220:                                     //'\'
        break;
    }
  });

  $('#baseCanvas').mousedown(function(event) {
    for (var i=0; i<scene.numPolyItems; i++) {
      if (scene.polyItems[i].pointIntersects([px, py])) {
        scene.dragItem = new dragItem(scene, [px, py], scene.polyItems[i].pIdx);
        break;
      }
    }
  });

  $('#baseCanvas').mousemove(function(event) {
    px = -canvas.width()/2 + event.clientX - offset.left;
    py = canvas.height()/2 - event.clientY + offset.top;
    if (scene.dragItem != null) {
      scene.dragItem.dragP = [px, py];
    }
  });

  $('#baseCanvas').mouseup(function(event) {
    scene.deleteDrag = true;
  });
});
