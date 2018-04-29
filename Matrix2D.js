function Matrix2D(
  n11, n12, n13, 
  n21, n22, n23, 
  n31, n32, n33 
){
  this.n11 = n11 || 1;
  this.n12 = n12 || 0;
  this.n13 = n13 || 0;
  this.n21 = n21 || 0;
  this.n22 = n22 || 1;
  this.n23 = n23 || 0;
  this.n31 = n31 || 0;
  this.n32 = n32 || 0;
  this.n33 = n33 || 1;
}

Matrix2D.prototype.clone = function()
{
  return new Matrix2D(this.n11, this.n12, this.n13, this.n21, this.n22, this.n23, this.n31, this.n32, this.n33);
}

Matrix2D.prototype.concat = function(m)
{
  var values = {};
  
  values.n11 = this.n11 * m.n11 + this.n12 * m.n21 + this.n13 * m.n31;
  values.n12 = this.n11 * m.n12 + this.n12 * m.n22 + this.n13 * m.n32;
  values.n13 = this.n11 * m.n13 + this.n12 * m.n23 + this.n13 * m.n33;
               
  values.n21 = this.n21 * m.n11 + this.n22 * m.n21 + this.n23 * m.n31;
  values.n22 = this.n21 * m.n12 + this.n22 * m.n22 + this.n23 * m.n32;
  values.n23 = this.n21 * m.n13 + this.n22 * m.n23 + this.n23 * m.n33;
               
  values.n31 = this.n31 * m.n11 + this.n32 * m.n21 + this.n33 * m.n31;
  values.n32 = this.n31 * m.n12 + this.n32 * m.n22 + this.n33 * m.n32;
  values.n33 = this.n31 * m.n13 + this.n32 * m.n23 + this.n33 * m.n33;
  
  this.initialize(values);
}

Matrix2D.prototype.initialize = function(values)
{
  for(var i in values) this[i] = values[i];
}

Matrix2D.prototype.createBox = function(scalex, scaley, scalez, rotationx, rotationy, rotationz, tx, ty, tz)
{
  this.identity();
  if (rotationx != 0) this.rotateX(rotationx);
  if (rotationy != 0) this.rotateY(rotationy);
  if (rotationz != 0) this.rotateZ(rotationz);
  if (scalex != 1 || scaley != 1 || scalez != 1) this.scale(scalex, scaley, scalez);
  if (tx != 0 || ty != 0 || tz != 0) this.translate(tx, ty, tz);
}

Matrix2D.prototype.identity = function()
{
  this.initialize({n11:1, n12:0, n13:0, n21:0, n22:1, n23:0, n31:0, n32:0, n33:1});
}

Matrix2D.prototype.rotate = function(angle)
{
  var sin = Math.sin(angle);
  var cos = Math.cos(angle);
  
  this.concat(new Matrix2D(
    cos, sin, 0,
    -sin, cos, 0,
    0, 0, 1)
  );
}

Matrix2D.prototype.scale = function(sx, sy)
{
  this.concat(new Matrix2D(
    sx, 0, 0, 
    0, sy, 0, 
    0, 0, 1)
  );
}

Matrix2D.prototype.translate = function(dx, dy)
{
  this.n31 += dx;
  this.n32 += dy;
}

Matrix2D.prototype.transformArray=function(arr)
{
    var rVal=[];

	var numPoints=arr.length/2;
	
	for(var i=0;i<numPoints;i++)
	{
		var i2=i*2;
		var x=arr[i2];
		var y=arr[i2+1];
		
		rVal[i2]=this.n11*x+this.n21*y+this.n31;
		rVal[i2+1]=this.n12*x+this.n22*y+this.n32;
	}
	
	return rVal;
}

Matrix2D.prototype.toString = function()
{
  return this.n11+","+this.n12+","+this.n13+","+this.n14+","+
    this.n21+","+this.n22+","+this.n23+","+this.n24+","+
    this.n31+","+this.n32+","+this.n33+","+this.n34+","+
    this.n41+","+this.n42+","+this.n43+","+this.n44;
}
