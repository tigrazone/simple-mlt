// Copyright (c) 2014 hole
// Copyright (c) 2015, 2016 tigra
// This software is released under the MIT License (http://kagamin.net/hole/license.txt).
// A part of this software is based on smallpt (http://www.kevinbeason.com/smallpt/) and
// released under the MIT License (http://kagamin.net/hole/smallpt-license.txt).
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <stack>
#include <ctime>
#include <cstring>

// OpenMP
#include <omp.h>

#define UINT_MAX ((unsigned int)-1)

const double PI = 3.14159265358979323846;
const double PI1_ = 1.0f/PI;
const double PI2 = PI + PI;
const double PI4 = PI2 + PI2;

const double INF = 1e20;
const double EPS = 1e-6;
const double MaxDepth = 5;

clock_t startT1,end1, end_t;


// *** その他の関数 ***
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; } 
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); } 

// Xorshift
// from Wikipedia
unsigned int xor128(void) { 
	static unsigned int x = 123456789;
	static unsigned int y = 362436069;
	static unsigned int z = 521288629;
	static unsigned int w = 88675123; 
	unsigned int t;

	t = x ^ (x << 11);
	x = y; y = z; z = w;
	return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)); 
}

const double _ui12345_ = 1.0/(UINT_MAX);

inline double rand01() { return (double)xor128() * _ui12345_; }

// *** データ構造 ***
struct Vec {
	double x, y, z;
	Vec(const double x_ = 0, const double y_ = 0, const double z_ = 0) : x(x_), y(y_), z(z_) {}
	inline Vec operator+(const Vec &b) const {return Vec(x + b.x, y + b.y, z + b.z);}
	inline Vec operator+(const double b) const {return Vec(x + b, y + b, z + b);}
	inline Vec operator-(const double b) const {return Vec(x - b, y - b, z - b);}
	inline Vec operator-(const Vec &b) const {return Vec(x - b.x, y - b.y, z - b.z);}
	inline Vec operator*(const double b) const {return Vec(x * b, y * b, z * b);}
	inline Vec operator/(const double b) const {const double b1 = 1.0 / b; return Vec(x * b1, y * b1, z * b1);}
	inline const double LengthSquared() const { return x*x + y*y + z*z; }
	inline const double Length() const { return sqrt(LengthSquared()); }
};
inline Vec operator*(double f, const Vec &v) { return v * f; }
inline Vec Normalize(const Vec &v) { return v / v.Length(); } //6* 1/ 1sqrt
// 要素ごとの積をとる
inline const Vec Multiply(const Vec &v1, const Vec &v2) {
	return Vec(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}
inline const double Dot(const Vec &v1, const Vec &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline const Vec Neg(const Vec &v1) {
	return Vec(-v1.x, -v1.y, -v1.z);
}

inline const Vec Cross(const Vec &v1, const Vec &v2) {
	return Vec((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
}
typedef Vec Color;
const Color BackgroundColor(0.0, 0.0, 0.0);


void save_hdr_file(const std::string &filename, const Color* image, const int width, const int height);
void print_time(clock_t endT,clock_t startT);

struct Ray {
	Vec org, dir;
	Ray(const Vec org_, const Vec &dir_) : org(org_), dir(dir_) {}
};

enum ReflectionType {
	DIFFUSE,    // 完全拡散面。いわゆるLambertian面。
	SPECULAR,   // 理想的な鏡面。
	REFRACTION, // 理想的なガラス的物質。
};

struct Sphere {
	double radius,r2;
	Vec position;
	Color emission, color;
	ReflectionType ref_type;
	Vec bbx_min,bbx_max;

	Sphere(const double radius_, const Vec &position_, const Color &emission_, const Color &color_, const ReflectionType ref_type_) :
	  radius(radius_), position(position_), emission(emission_), color(color_), ref_type(ref_type_) {
	  r2=radius_*radius_;
	  bbx_min=position_ - radius_;
	  bbx_max=position_ + radius_;
	  }
	// 入力のrayに対する交差点までの距離を返す。交差しなかったら0を返す。
	const double intersect(const Ray &ray) {
	
	/*
with 65s
without 55s	
*/
/*
	//x	
	if (ray.dir.x<0.0)
	{
		if (ray.org.x<bbx_min.x)
			return 0.0;
	}
	else
	{
		if (ray.org.x>bbx_max.x)
			return 0.0;
	}
	
	//y	
	if (ray.dir.y<0.0)
	{
		if (ray.org.y<bbx_min.y)
			return 0.0;
	}
	else
	{
		if (ray.org.y>bbx_max.y)
			return 0.0;
	}
	
	//z	
	if (ray.dir.z<0.0)
	{
		if (ray.org.z<bbx_min.z)
			return 0.0;
	}
	else
	{
		if (ray.org.z>bbx_max.z)
			return 0.0;
	}
	 */
		Vec o_p = position - ray.org;
		
		const double b = Dot(o_p, ray.dir), det = b * b - Dot(o_p, o_p) + r2;
		//7*
		if (det >= 0.0) {
			const double sqrt_det = sqrt(det);
			const double t1 = b - sqrt_det, t2 = b + sqrt_det;
			if (t1 > EPS)		return t1;
			else if(t2 > EPS)	return t2;
		}
		return 0.0;
		
	}
};

// *** レンダリングするシーンデータ ****
// from small ppt
//9 spheres
//radius, position, emission, color, mat

Sphere spheres[] = {
	Sphere(5.0, Vec(50.0, 75.0, 81.6),Color(12,12,12), Color(), DIFFUSE),//照明
	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Color(), Color(0.75, 0.25, 0.25),DIFFUSE),// 左
	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Color(), Color(0.25, 0.25, 0.75),DIFFUSE),// 右
	Sphere(1e5, Vec(50,40.8, 1e5), Color(), Color(0.75, 0.75, 0.75),DIFFUSE),// 奥
	Sphere(1e5, Vec(50,40.8,-1e5+170), Color(), Color(), DIFFUSE),// 手前
	Sphere(1e5, Vec(50, 1e5, 81.6), Color(), Color(0.75, 0.75, 0.75),DIFFUSE),// 床
	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Color(), Color(0.75, 0.75, 0.75),DIFFUSE),// 天井
	Sphere(16.5,Vec(27,16.5,47), Color(), Color(1,1,1)*.99, SPECULAR),// 鏡
	Sphere(16.5,Vec(73,16.5,78), Color(), Color(1,1,1)*.99, REFRACTION),//ガラス
};
const int LightID = 0;

/*
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFFUSE),//Left
  Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFFUSE),//Rght
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFFUSE),//Back
  Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFFUSE),//Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFFUSE),//Botm
  Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFFUSE),//Top
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPECULAR),//Mirr
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFRACTION),//Glas
  Sphere(1.5, Vec(50,81.6-16.5,81.6),Vec(4,4,4)*100,  Vec(), DIFFUSE),//Lite
};
const int LightID = 8;
*/

// *** レンダリング用関数 ***
// シーンとの交差判定関数
inline bool intersect_scene(const Ray &ray, double *t, int *id) {
	const double n = sizeof(spheres) / sizeof(Sphere);
	*t  = INF;
	*id = -1;
	for (int i = 0; i < int(n); i ++) {
		double d = spheres[i].intersect(ray);
		if (d > 0.0 && d < *t) {
			*t  = d;
			*id = i;
		}
	}
	return *t < INF;
}

// A Simple and Robust Mutation Strategy for the Metropolisを参照。
// Kelemen style MLT用データ構造
// Kelemen styleではパス生成に使う乱数の空間で変異させたりする。
// その一つ一つのサンプルのデータ構造。
struct PrimarySample {
	int modify_time;
	double value;
	PrimarySample() {
		modify_time = 0;
		value = rand01();
	}
};

const double s2 = 1.0 / 16.0;

double s1 = 1.0 / 256.0; 
double s12 = s1 / s2;
double s112 = s1 / (s12 + 1.0);

// Kelemen MLTにおいて、パス生成に使う各種乱数はprimary spaceからもってくる。
// PrimarySample()を通常のrand01()の代わりに使ってパス生成する。今回は普通のパストレースを使った。（双方向パストレ等も使える）
// Metropolis法なので現在の状態から次の状態へと遷移をするがこの遷移の空間が従来のワールド空間ではなく
// パスを生成するのに使われる乱数の空間になっている。
// 乱数空間における変異（Mutate()）の結果を使って再び同じようにパストレでパスを生成するとそのパスは自然に変異後のパスになっている。
struct KelemenMLT {
private:	
		
	// LuxR}erから拝借してきた変異関数
	inline double Mutate(const double  x)
	{
	/*
		const double r = rand01();
		const double s1 = 1.0 / 512.0, s2 = 1.0 / 16.0;
		const double dx = s1 / (s1 / s2 + fabsf(2.0 * r - 1.0)) - s1 / (s1 / s2 + 1.0);
		// 6/
		if (r < 0.5) {
			const double x1 = x + dx;
			return (x1 < 1.0) ? x1 : x1 - 1.0;
		} else {
			const double x1 = x - dx;
			return (x1 < 0.0) ? x1 + 1.f : x1;
		}		
		*/
		
		const double r = rand01();
		const double dx = s1 / (s12 + fabsf(r + r - 1.0)) 
		- s112;
		//tigra: little speedup  1/  -5/  15% faster
		
		if (r < 0.5) {
			const double x1 = x + dx;
			return (x1 < 1.0) ? x1 : x1 - 1.0;
		} else {
			const double x1 = x - dx;
			return (x1 < 0.0) ? x1 + 1.f : x1;
		}
	}
public:
	// 論文で使われているものとほぼ同じ
	int global_time;
	unsigned char large_step;
	int large_step_time;
	int used_rand_coords;

	std::vector<PrimarySample> primary_samples;
	std::stack<PrimarySample> primary_samples_stack;

	KelemenMLT() {
		global_time = large_step = large_step_time = used_rand_coords = 0;
		primary_samples.resize(128);
	}
	void InitUsedRandCoords() {
		used_rand_coords = 0;
	}

	// 通常の乱数のかわりにこれを使う
	// 論文にのっているコードとほぼ同じ
	// Используйте его вместо обычного случайного числа
	// Код, который ехал в статье о том же
	inline double NextSample() {
		if (primary_samples.size() <= used_rand_coords) {
			primary_samples.resize(primary_samples.size() * 1.5); // 拡張する расширить
		}

		if (primary_samples[used_rand_coords].modify_time < global_time) {
			//if (large_step > 0) {
			if (large_step) {
				primary_samples_stack.push(primary_samples[used_rand_coords]);
				primary_samples[used_rand_coords].modify_time = global_time;
				primary_samples[used_rand_coords].value = rand01();
			} else {
				if (primary_samples[used_rand_coords].modify_time < large_step_time) {
					primary_samples[used_rand_coords].modify_time = large_step_time;
					primary_samples[used_rand_coords].value = rand01();
				}

				while (primary_samples[used_rand_coords].modify_time < global_time - 1) {
					primary_samples[used_rand_coords].value = Mutate(primary_samples[used_rand_coords].value);
					primary_samples[used_rand_coords].modify_time ++;
				}
				
				primary_samples_stack.push(primary_samples[used_rand_coords]);
				primary_samples[used_rand_coords].value = Mutate(primary_samples[used_rand_coords].value);
				primary_samples[used_rand_coords].modify_time = global_time;
			}
		}

		used_rand_coords ++;
		return primary_samples[used_rand_coords - 1].value;
	}
};


double luminance(const Color &color) {
	return Dot(Vec(0.2126, 0.7152, 0.0722), color);
}

// 光源上の点をサンプリングして直接光を計算する
Color direct_radiance_sample(const Vec &v0, const Vec &normal, const int id, 
KelemenMLT &mlt) 
{
	// 光源上の一点をサンプリングする
	const double r1 = PI2 * mlt.NextSample();
	const double r2 = 1.0 - 2 * mlt.NextSample();
	
	//40s vs 48s 3m
	const double r22 = sqrt(1.0 - r2*r2);
	
	//tigra: Light tracing. в реальном приложении надо выбирать с помощью техники из psychopath
	const Vec light_pos = spheres[LightID].position + 
				((spheres[LightID].radius + EPS) * 
				Vec(r22 * cos(r1), 
				r22 * sin(r1), 
				r2));
	
	// サンプリングした点から計算
	const Vec light_normal = Normalize(light_pos - spheres[LightID].position);
	const Vec light_dir = 
	Normalize
	(light_pos - v0);
	
	const double dot0 = Dot(normal, light_dir);
	const double dot1 = -Dot(light_normal, light_dir);

	if (dot0 >= 0 && dot1 >= 0) 
	{
		
		double t; // レイからシーンの交差 位置までの距離
		int id_; // 交差したシーン内オブジェクトのID
		//light cache here add
		if(!intersect_scene(Ray(v0, light_dir), &t, &id_))			
			return Color();
		
	//tigra
	const double dist2 = (light_pos - v0).LengthSquared();
	
		if (fabs(sqrt(dist2) - t) < 1e-3) {	
		
	
		const double G = 4.0f *dot0 * dot1 / dist2;
		
			return Multiply(spheres[id].color, spheres[LightID].emission) *
					( G * ( ( 
					//pow(spheres[LightID].radius, 2.0)
					//tigra
					spheres[LightID].r2
					)
					));
		}
	}
	return Color();
}

// ray方向からの放射輝度を求める
// ただし、rand01()の代わりにKelemenMLT::NextSample()を使う。
// それ以外は普通のパストレースと同じ。
// ただし、各点で明示的に光源の影響をサンプリングして求める。（そのほうがMLTとの相性がよい？）
// あるいは双方向パストレーシングにするのもよいと思う。
Color radiance(const Ray &ray, const int depth, KelemenMLT &mlt) {
	double t; // レイからシーンの交差位置までの距離
	int id;   // 交差したシーン内オブジェクトのID
	if (!intersect_scene(ray, &t, &id))
		return BackgroundColor;

	const Sphere &obj = spheres[id];
	const Vec hitpoint = ray.org + t * ray.dir; // 交差位置
	const Vec normal  = Normalize(hitpoint - obj.position); // 交差位置の法線
								
	// 色の反射率最大のものを得る。ロシアンルーレットで使う。
	// ロシアンルーレットの閾値は任意だが色の反射率等を使うとより良い。
	double russian_roulette_probability = std::max(obj.color.x, 
										  std::max(obj.color.y, 
												   obj.color.z)
												  );
	// 一定以上レイを追跡したらロシアンルーレットを実行し追跡を打ち切るかどうかを判断する
	if (depth > MaxDepth) {
		if (mlt.NextSample() >= russian_roulette_probability)
			return Color();
	} else
		russian_roulette_probability = 1.0; // ロシアンルーレット実行しなかった
	
	const double rr1 = 1.0f/russian_roulette_probability;

	switch (obj.ref_type) {
	case DIFFUSE: {
		// 直接光のサンプリングを行う
		// Я делать выборку из прямого света
		if (id != LightID) {
			
	const Vec orienting_normal = Dot(normal, ray.dir) < 0.0 ? normal : 
								Neg( normal); 
								// 交差位置の法線（物体からのレイの入出を考慮）
								
			const int shadow_ray = 4;
			const double sr1 = 1.0 / shadow_ray;
			Vec direct_light;
			for (int i = 0; i < shadow_ray; i ++) {
				direct_light = direct_light + 
					direct_radiance_sample(hitpoint, orienting_normal, id, mlt) 
					// / shadow_ray;
					* sr1;
			}

			// orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
			Vec w, u, v;
			
			w = orienting_normal;
			if (fabs(w.x) > 0.1)
			{
				//u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
				
				const double sq = 1.0 / sqrt (w.z * w.z + w.x * w.x);
				//u = Normalize(Vec(w.z, 0.0, -w.x)); //-6*
				u = Vec(w.z*sq, 0.0, -w.x*sq); //was 6* 1/ 1sqrt -> 4* 1/ 1sqrt  => -2*				
			
				v = Vec(w.y*u.z, w.z*u.x - w.x*u.z, -w.y*u.x);
				
			}
			else
			{
				//u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
				
				const double sq = 1.0 / sqrt (w.z * w.z + w.y * w.y);
				//u = Normalize(Vec(0.0, -w.z, w.y)); //-6*
				
				u = Vec(0.0, -w.z*sq, w.y*sq); //-2*
				v = Vec(w.y*u.z - w.z*u.y, -w.x*u.z, w.x*u.y);
			}
			
			//v = Cross(w, u); 
			
			//6* 6* 1/ 1sqrt 6* 6* 1/ 1sqrt => 24* 2/ 2sqrt =>   8* 1/ 1sqrt => -16* 1/ 1sqrt
			
			// コサイン項を使った重点的サンプリング
			// Средоточие выборки с использованием терминов косинуса
			const double r1 = PI2 * mlt.NextSample();
			const double r2 = mlt.NextSample(), r2s = sqrt(r2);
			Vec dir = 
			//tigra: not needed, speedup
			//Normalize
								((u * cos(r1) * r2s + 
								 v * sin(r1) * r2s + 
								 w * sqrt(1.0 - r2)
								));

			return (direct_light + Multiply(obj.color, 
					radiance(Ray(hitpoint, dir), depth+1, mlt))) * rr1 
					;
		} else if (depth == 0) {
			return obj.emission;
		} else
			return Color();
	} break;
	case SPECULAR: {
		// 完全鏡面なのでレイの反射方向は決定的。
		// ロシアンルーレットの確率で除算するのは上と同じ。
		// С полным направлении зеркального отражения Ray решающей.
		// То же самое и с делится на вероятность русскую рулетку.
		double lt;
		int lid;
		
		const double dnd = Dot(normal, ray.dir);
		
		Ray reflection_ray = Ray(hitpoint, ray.dir - normal * (dnd + dnd) );
		intersect_scene(reflection_ray, &lt, &lid);
		Vec direct_light;
		if (lid == LightID)
			direct_light = spheres[LightID].emission;

		return (direct_light + Multiply(obj.color, radiance(reflection_ray, depth+1, mlt))) * rr1;
	} break;
	case REFRACTION: {
		const double dnd = Dot(normal, ray.dir);
		
		Ray reflection_ray = Ray(hitpoint, ray.dir - normal * ( dnd + dnd));
		
		// 反射方向からの直接光
		// Прямой свет от направления отражения
		double lt;
		int lid;
		intersect_scene(reflection_ray, &lt, &lid);
		Vec direct_light;
		if (lid == LightID)
			direct_light = spheres[LightID].emission;
		
	const Vec orienting_normal = dnd < 0.0 ? normal : 
								Neg( normal); 
								// 交差位置の法線（物体からのレイの入出を考慮）

		bool into = Dot(normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか 
					// Если луч выходит из объекта или от входа

		// Snellの法則
		const double nc = 1.0; // 真空の屈折率
		const double nt = 1.5; // オブジェクトの屈折率
		const double nnt = into ? nc / nt : nt / nc;
		const double ddn = Dot(ray.dir, orienting_normal);
		const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
		
		if (cos2t < 0.0) { // 全反射した
			return (direct_light + Multiply(obj.color, 
			(radiance(reflection_ray, depth+1, mlt)))) * rr1 
			;
		}
		
		
		//tigra: 01 little bit optimize - remove *1.0
		double aa =(ddn * nnt + sqrt(cos2t));
		if(!into) aa=-aa;
		
		// 屈折していく方向
		Vec tdir = 
		Normalize
		(ray.dir * nnt - normal * aa);
		
		/*
		// 屈折していく方向
		Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * 
		(ddn * nnt + sqrt(cos2t)));
		*/
		
		// SchlickによるFresnelの反射係数の近似
		const double a = nt - nc, b = nt + nc;
		const double R0 = (a * a) / (b * b);
		const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));	
		
		//const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
		
		//tigra: 02 pow(c, 5.0) -> 3*
		//const double cc = c*c;
		const double Re = R0 + (1.0 - R0) * c * c * c * c * c;
		
		//tigra: 01+02  3% speedup
		
		const double Tr = 1.0 - Re; // 屈折光の運ぶ光の量
		const double probability  = 0.25 + 0.5 * Re;
		const double probability1  = 1.0 - probability;
		
		// 屈折方向からの直接光
		Ray refraction_ray = Ray(hitpoint, tdir);
		intersect_scene(refraction_ray, &lt, &lid);
		
		Vec direct_light_refraction;
		
		//tigra: we HIT the LIGHT!!!!!!
		if (lid == LightID)
			direct_light_refraction = spheres[LightID].emission;

		// 一定以上レイを追跡したら屈折と反射のどちらか一方を追跡する。（さもないと指数的にレイが増える）
		// ロシアンルーレットで決定する。
		// Для отслеживания либо отражение и преломление После отслеживания более определенного Луча. (В противном случае экспоненциальных нуля увеличивается)
		// Я буду определяется русскую рулетку.
		if (depth > 2) {
			if (mlt.NextSample() < probability) { // 反射
			const double prbblty=probability;
			
				return  Multiply(obj.color, (direct_light + radiance(reflection_ray, depth+1, mlt)) * Re)
					/ (
					prbblty
					* rr1
					);
			} else { // 屈折
			const double prbblty=probability1;
			
				return  Multiply(obj.color, (direct_light_refraction + radiance(refraction_ray, depth+1, mlt)) * Tr)
					/ (
					prbblty
					* rr1
					);
			}
		} else { // 屈折と反射の両方を追跡
			return Multiply(obj.color, (direct_light + radiance(reflection_ray, depth+1, mlt)) * Re
				                  + (direct_light_refraction + radiance(refraction_ray, depth+1, mlt)) * Tr) * rr1;
		}
	} break;
	}

	return Color();
}

// 上のパストレで生成したパスを保存しておく
// Чтобы сохранить тот путь, который вы сгенерированный выше
struct PathSample {
	int x, y;

	Color F;
	double weight;
	PathSample(const int x_ = 0, const int y_ = 0, const Color &F_ = Color(), const double weight_ = 1.0) :
	x(x_), y(y_), F(F_), weight(weight_) {}
};

#define R1(x) ((x)+(x))

// MLTのために新しいパスをサンプリングする関数。 Функция сэмплирования новый путь для
// 今回はradiance()（パストレ）を使ったが何でもいい。
PathSample generate_new_path(const Ray &camera, const Vec &cx, const Vec &cy, const int width, const int height, const double &width1, const double &height1, KelemenMLT &mlt, int x, int y) {
	double weight = 4.0;
	
	//-1 всегда
	if (x < 0) {
		weight *= width;
		x = mlt.NextSample() * width;
		if (x == width)
			x = 0;
	}
	
	//-1 всегда
	if (y < 0) {
		weight *= height;
		y = mlt.NextSample() * height;
		if (y == height)
			y = 0;
	}
	
	int sx = mlt.NextSample() < 0.5 ? 0 : 1;
	int sy = mlt.NextSample() < 0.5 ? 0 : 1;
	
	// テントフィルターによってサンプリング
	// ピクセル範囲で一様にサンプリングするのではなく、ピクセル中央付近にサンプルがたくさん集まるように偏りを生じさせる
	
	// Оцифровывается палаточном фильтра
	// Вместо равномерно пробы в диапазоне пикселей, чтобы производить смещение в сбора образца А много в непосредственной близости от центра пикселя
	//const double r01 = mlt.NextSample();
	//const double r1 = r01 + r01;
	
	const double r1 = 2 * mlt.NextSample();
	const double dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
	//const double dx = R1(r01) < 1.0 ? sqrt(R1(r01)) - 1.0 : 1.0 - sqrt(2.0 - R1(r01));
	//const double dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
				 
	const double r2 = 2 * mlt.NextSample(),
				 dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
				 
				 /*
	const double r02 = mlt.NextSample();
	//const double r2 = r02 + r02;
	
	const double 
				 dy = R1(r02) < 1.0 ? sqrt(R1(r02)) - 1.0 : 1.0 - sqrt(2.0 - R1(r02));
				 //dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
				 */
	//     / 2.0 -> *0.5
	/*			 
	Vec dir = cx * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5) +
				cy * (((sy + 0.5 + dy) / 2.0 + y) / height- 0.5) + camera.dir;
	*/
	
	Vec dir = cx * (((sx + 0.5 + dx) *0.5 + x) * width1 - 0.5) +
	//Vec dir = cx * (((sx + 0.5 + dx) *0.5 + x) / width - 0.5) +
				cy * (((sy + 0.5 + dy) *0.5 + y) * height1- 0.5) + camera.dir;
	//			cy * (((sy + 0.5 + dy) *0.5 + y) / height- 0.5) + camera.dir;
	
	const Ray ray = Ray(camera.org + dir * 130.0, 
	Normalize
	(dir));

	Color c = radiance(ray, 0, mlt);

	return PathSample(x, y, c,  weight);
}

// MLTする
// 画像平面上で大域的に変異させる глобально мутация на плоскости изображения
void render_mlt(const int mlt_num, const long mutation, Color *image, const Ray &camera, const Vec &cx, const Vec &cy, const int width, const int height) {
	// MLTを別々に並列に走らせて結果をマージする Чтобы объединить результат давали отдельно параллельно

	const long wh0 = width * height;		
	const double mlt_num1 = 1.0f / mlt_num;
	
	const double width1 = 1.0f / width;
	const double height1 = 1.0f / height;
	
		// 1/16 = 0.0625 = 6.25%
		// 1/8 = 0.125 = 12.5%
	
		//Получить первый путь будет использоваться с этого пути в MLT.
		//const long SeedPathMax = ((wh0 / 16)+(wh0 / 8)) / 2; // 適当に多めの数 //Количество соответствующим много
		// среднее 9.375%
		
		//const long SeedPathMax = ((wh0 >> 4)+(wh0 >> 3) + 16 ) >> 1; // 適当に多めの数 //Количество соответствующим много
		
		// /16  6.25%
		//const long SeedPathMax = (wh0 >> 4); // 適当に多めの数 //Количество соответствующим много
		const long SeedPathMax = wh0 / 10; // 適当に多めの数 //Количество соответствующим много
		//const long SeedPathMax = mutation >> 4; // /3=8 4=16 5=32 6=64 7=128 8=256   /32 = 3%   /16 = 6.6%   /8 = 12.5%
		//int SeedPathMax = mutation ; // 適当に多めの数 //Количество соответствующим много
		
		/*
		SeedPathMax*=0.05; // 5% новых
		//SeedPathMax*=2; // 10% новых
		SeedPathMax+=SeedPathMax; // 10% новых
		*/
		
		
	int num_threads = 0;
	int num_cores = omp_get_num_procs();

	if(num_threads<1)
		num_threads  = std::max(1, num_cores )*4;
	
	printf("num_threads=%d num_cores=%d\n", num_threads, num_cores);
		
    // Set number of used threads
    omp_set_num_threads(num_threads);
	
	// OpenMP
	 // omp_lock_t lock0;
	 // omp_init_lock(&lock0);
	 
 #pragma omp parallel for schedule(dynamic)
	for (int mi = 0; mi < mlt_num; mi ++) {
		if(mlt_num>1)
			printf("\nMLT pass: %d\n",mi+1);
	
		std::vector<Color> tmp_image;
		tmp_image.resize(wh0);

		KelemenMLT mlt;
		// たくさんパスを生成する。
		// このパスからMLTで使う最初のパスを得る。(Markov Chain Monte Carloであった）
		//генерировать путь много.
		/*
		if (SeedPathMax <= 0)
			SeedPathMax = 1;
		*/
			
		printf("creating SeedPaths %d...\n",SeedPathMax);
		
		end_t=clock();

		std::vector<PathSample> seed_paths(SeedPathMax);
		double sumI = 0.0;
		mlt.large_step = 1;
		for (long i = 0; i < SeedPathMax; i ++) {
			mlt.InitUsedRandCoords();
			PathSample sample = generate_new_path(camera, cx, cy, width, height, width1, height1, mlt, -1, -1);
			mlt.global_time ++;
			while (!mlt.primary_samples_stack.empty()) // スタック空にする //стек пустой
				mlt.primary_samples_stack.pop();

			sumI += luminance(sample.F);
			seed_paths[i] = sample;
		}
		end1=clock();
				
				print_time(end1, end_t);
		
		
			
		// 最初のパスを求める。輝度値に基づく重点サンプリングによって選んでいる。
		//Мы просим первый путь. Вы выбрали с помощью выборки по значимости на основе значения яркости.
		const double rnd = rand01() * sumI;
		int selecetd_path = 0;
		double accumulated_importance = 0.0;
		for (long i = 0; i < SeedPathMax; i ++) {
			accumulated_importance += luminance(seed_paths[i].F);
			if (accumulated_importance >= rnd) {
				selecetd_path = i;
				break;
			}
		}
	
		// 論文参照
		//const double b = sumI / SeedPathMax;
		
			const double b1 = SeedPathMax / sumI ;
		
		const double p_large = 0.5;
		const long M = mutation;
		
		long idx;
		
		const double M1 = 1.0f / M;
		double olda_time = clock();
		
		//const long M2 = M / 100; //1%
		//compiler do it faster!
		
		int accept = 0, reject = 0;	
		PathSample old_path = seed_paths[selecetd_path];
		int progress = 0;
		
		for (long i = 0; i < M ; i ++) {
			if ((i + 1) % (M / 100) == 0 || progress==99) {
				end1=clock();
					
				progress ++; //1%
				
				//printf("delta_time=%.2f\n",end1-olda_time);
				
				if(end1-olda_time>3000 || progress==100)
				{		
					olda_time=clock();
				std::cout << progress << "%  " ;
				std::cout << "Accept: " << accept << " Reject: " << reject << " Rate: " << (100.0 * accept / (accept + reject)) << "%" << std::endl;
				
				
				print_time(end1,startT1);
				}
			}

			// この辺も全部論文と同じ（Next()）
			//Эта сторона также же, как и всей работы
			mlt.large_step = rand01() < p_large;
			mlt.InitUsedRandCoords();	
			
			PathSample new_path = generate_new_path(camera, cx, cy, width, height, width1, height1, mlt, -1, -1);
			
			/*
			double a = std::min(1.0, luminance(new_path.F) / luminance(old_path.F));
			const double new_path_weight = (a + mlt.large_step) / (luminance(new_path.F) / b + p_large) / M;
			const double old_path_weight = (1.0 - a) / (luminance(old_path.F) / b + p_large) / M;
			
			7/
			1/ вне цикла
			*/
			
			//const double lum_n = luminance(new_path.F);
			const double lum_n = luminance(new_path.F)*b1;
			//const double lum_o = luminance(old_path.F);
			const double lum_o = luminance(old_path.F)*b1;
						
			const double a = std::min(1.0, lum_n / lum_o);
			const double M1a = M1 * a;
			
			/*
			double new_path_weight;
			
			if(mlt.large_step)
				new_path_weight = (M1a + M1 * mlt.large_step) / (lum_n + p_large) ;
			else
				new_path_weight = (M1a) / (lum_n + p_large) ;
			*/
			
			const double new_path_weight = (M1a + mlt.large_step * M1 ) / (lum_n + p_large) ;
			//const double new_path_weight = (M1a + ((mlt.large_step)?M1 :0) ) / (lum_n + p_large) ;
			//const double new_path_weight = (M1a + M1 * mlt.large_step) / (lum_n * b1 + p_large) ;
			//const double new_path_weight = (a + mlt.large_step) / (lum_n * b1 + p_large) * M1;
			
			const double old_path_weight = (M1 - M1a) / (lum_o + p_large) ;
			//const double old_path_weight = (M1 - M1a) / (lum_o * b1 + p_large) ;
			//const double old_path_weight = (1.0 - a) / (lum_o * b1 + p_large) * M1;
			/*
			3/ 3*
			1/ вне цикла
			*/
			
			idx = new_path.y * width + new_path.x;			
			tmp_image[idx] = tmp_image[idx] + (new_path.weight * new_path_weight) * new_path.F;
			
			idx = old_path.y * width + old_path.x;			
			tmp_image[idx] = tmp_image[idx] + (old_path.weight * old_path_weight) * old_path.F;
				
				//a*=2; //tigra: this make 60% accept vs ~50%
				//a+=a; //tigra: this make 60% accept vs ~50%
				
			//if (rand01() < (a+a)) 
			if (rand01() < (a*1.5)) 
			{ 
			// 受理 //принятие
				accept ++;
				old_path = new_path;
				if (mlt.large_step)
					mlt.large_step_time = mlt.global_time;
				mlt.global_time ++;
				while (!mlt.primary_samples_stack.empty()) // スタック空にする //стек пустой
					mlt.primary_samples_stack.pop();
			} else { // 棄却 //отказ
				reject ++;
				idx = mlt.used_rand_coords - 1;
				while (!mlt.primary_samples_stack.empty()) {
					mlt.primary_samples[idx --] = mlt.primary_samples_stack.top();
					mlt.primary_samples_stack.pop();
				}
			}
		}
		
		// OpenMP
		 // omp_set_lock(&lock0);
		 
		//for(long i = 0; i < width * height; i ++) {
		for(long i = 0; i < wh0; i ++) {
			//image[i] = image[i] + tmp_image[i] / mlt_num;
			image[i] = image[i] + tmp_image[i] * mlt_num1;
		}
		 // omp_unset_lock(&lock0);
	}
	// OpenMP
	 // omp_destroy_lock(&lock0);
}

// *** .hdrフォーマットで出力するための関数 ***
struct HDRPixel {
	unsigned char r, g, b, e;
	HDRPixel(const unsigned char r_ = 0, const unsigned char g_ = 0, const unsigned char b_ = 0, const unsigned char e_ = 0) :
	r(r_), g(g_), b(b_), e(e_) {};
	unsigned char get(int idx) {
		switch (idx) {
		case 0: return r;
		case 1: return g;
		case 2: return b;
		case 3: return e;
		} return 0;
	}

};

// doubleのRGB要素を.hdrフォーマット用に変換
HDRPixel get_hdr_pixel(const Color &color) {
	double d = std::max(color.x, std::max(color.y, color.z));
	if (d <= 1e-32)
		return HDRPixel();
	int e;
	double m = frexp(d, &e); // d = m * 2^e
	d = m * 256.0 / d;
	return HDRPixel(color.x * d, color.y * d, color.z * d, e + 128);
}

// 書き出し用関数
void save_hdr_file(const std::string &filename, const Color* image, const int width, const int height) {
	FILE *fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		std::cerr << "Error: " << filename << std::endl;
		return;
	}
	// .hdrフォーマットに従ってデータを書きだす
	// ヘッダ
	unsigned char ret = 0x0a;
	fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
	fprintf(fp, "# Made with 100%% pure HDR Shop%c", ret);
	fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
	fprintf(fp, "EXPOSURE=1.0000000000000%c%c", ret, ret);

	// 輝度値書き出し
	fprintf(fp, "-Y %d +X %d%c", height, width, ret);
	for (int i = height - 1; i >= 0; i --) {
		std::vector<HDRPixel> line;
		for (int j = 0; j < width; j ++) {
			HDRPixel p = get_hdr_pixel(image[j + i * width]);
			line.push_back(p);
		}
		fprintf(fp, "%c%c", 0x02, 0x02);
		fprintf(fp, "%c%c", (width >> 8) & 0xFF, width & 0xFF);
		for (int i = 0; i < 4; i ++) {
			for (int cursor = 0; cursor < width;) {
				const int cursor_move = std::min(127, width - cursor);
				fprintf(fp, "%c", cursor_move);
				for (int j = cursor;  j < cursor + cursor_move; j ++)
					fprintf(fp, "%c", line[j].get(i));
				cursor += cursor_move;
			}
		}
	}

	fclose(fp);
}

//tigra: time functions
void get_time_str(char * s,clock_t endT,clock_t startT)
{
	int hours,mins,secs;
		
		float elapsedTime = float((endT - startT) / CLOCKS_PER_SEC);
		
		hours=(int)floor(elapsedTime/3600);
		mins=(int)floor((elapsedTime-hours*3600)/60);
		secs=(int)floor(elapsedTime-hours*3600-mins*60);			

		sprintf(s,"%0.2fs",elapsedTime);		
		
			if(elapsedTime>=60)
				{
				 strcat(s,"(");
				 if(hours>0) sprintf(s,"%s%dh",s,hours);
				 if(mins>0) 
					{
					if(hours>0) strcat(s,":");
					sprintf(s,"%s%dmin:",s,mins);
					}
				 sprintf(s,"%s%02ds)",s,secs);
				}
}

void print_time(clock_t endT,clock_t startT)
{	
		char s[1024];

		get_time_str(s,endT,startT);
		printf("%s\n",s);
}

//tigra: writes k M G for number
char get_maga_giga_str(double &num)
{			
		double kilo=num/1000;
		double mega=kilo/1000;
		double giga=mega/1000;

		char c=' ';

		if(kilo>0.999f)
		{
			c='k';
			num=kilo;
			if(mega>0.999f)
			{
				if(giga>0.999f)
				{
					c='G';
					num=giga;
				}
				else
					{
					c='M';
					num=mega;
				}

			}
		}

  return c;
}

void PrintHelp()
{
	printf("-mutations muts\tmuts mutations\n");
	printf("-mpp mpps\tmutations per pixel\n");
	printf("-passes passes\tnumber of MLT passes\n");
	printf("-width w\twidth of outout image w pixels\n");
	printf("-height h\theight of outout image h pixels\n");
}

//g++ -O3 -fopenmp simplemlt.cpp
int main(int argc, char **argv) {
	
	int itmp;
	int m_chng=0;
	long mutation;
	
	printf("SimpleMLT\nby Hole http://kagamin.net/hole and tigra http://thrlite.blogspot.com\n");
	printf("18 jun 2016\n");
	
	printf("-h\thelp\n");
	printf("/?\t\n");
	printf("--help\t\n");
	
	printf("\n");
	
	int width = 1024;
	int height = 768;
	//int mutation = 32 * width * height;
	
	long muts=16;
	
	int mlt_num = 1; // 何回MLTするか Сколько раз MLT
	
	//mutation=5*1000*1000;
	
	//tigra: read from command line
	if(argc>1)
	{
		// Load arguments
		for(int i=1; i<argc; i++)
		{
			std::string arg(argv[i]);

			// print help string (at any position)
			if(arg == "-h" || arg == "--help" || arg == "/?" || arg == "-?")
			{
				PrintHelp();
				return 0;
			}

			if(arg[0] != '-') // all our commands start with -
			{
				continue;
			}
			else if(arg == "-mpp")
			{
				if(++i == argc)
				{
					printf("Missing <mutations> argument, use default\n");
				}

				itmp=atoi(argv[i]);

				if(itmp < 1)
				{
					printf("Invalid <mutations> argument, use default\n");
				}
				else
					{
						muts=itmp;
					}
			}
			else if(arg == "-passes")
			{
				if(++i == argc)
				{
					printf("Missing <passes> argument, use default\n");
				}

				itmp=atoi(argv[i]);

				if(itmp < 1)
				{
					printf("Invalid <passes> argument, use default\n");
				}
				else
					{
						mlt_num=itmp;						
					}
			}
			else if(arg == "-mutations")
			{
				if(++i == argc)
				{
					printf("Missing <mutations> argument, use default\n");
				}

				itmp=atoi(argv[i]);

				if(itmp < 1)
				{
					printf("Invalid <mutations> argument, use default\n");
				}
				else
					{
						mutation=itmp;
						m_chng=1;
					}
			}
			else if(arg == "-width")
			{
				if(++i == argc)
				{
					printf("Missing <width> argument, use default\n");
				}
				
				itmp=atoi(argv[i]);

				if(itmp < 1)
				{
					printf("Invalid <width> argument, use default\n");
				}
				else
					width=itmp;
			}
			else if(arg == "-height")
			{
				if(++i == argc)
				{
					printf("Missing <height> argument, use default\n");
				}
				
				itmp=atoi(argv[i]);

				if(itmp < 1)
				{
					printf("Invalid <height> argument, use default\n");
				}
				else
					height=itmp;
			}
		}	
	}
	/*	
	else
	{
				PrintHelp();
				return 0;
			}
	*/
	
	if(!m_chng)
		{
			printf("Mutations per pass: %d\n",muts);
			
			mutation = muts * width * height;
		}
	
	if(mlt_num>1)
		{
			printf("MLT passes: %d\n",mlt_num);
		}
		
	int nnn = (int) sizeof(spheres) / sizeof(Sphere);	
	printf("objects: %d\n", nnn);
	
	double mmm=mutation;
	
	char kmg=get_maga_giga_str(mmm);
		
	//printf("muts=%dx%dx%d\n",muts,width,height);
	
	if(mlt_num>1)
	printf("mutations per pass %.2f%c\n", mmm,kmg);
	else
	printf("mutations %.2f%c\n", mmm,kmg);
	
	
	s1 = 1.0f / (width>>1); 
	s12 = s1 / s2;
	s112 = s1 / (s12 + 1.0f);
	

	// カメラ位置
	// положение камеры
	Ray camera(Vec(50.0, 52.0, 295.6), Normalize(Vec(0.0, -0.042612, -1.0)));
	// シーン内でのスクリーンのx,y方向のベクトル
	// Экран х на сцене, у направление вектора
	Vec cx = Vec(width * 0.5135 / height);
	Vec cy = Normalize(Cross(cx, camera.dir)) * 0.5135;
	Color *image = new Color[width * height];
	
	
	startT1=clock();
	
	render_mlt(mlt_num, mutation, image, camera, cx, cy, width, height);
	
	end1=clock();
	
	// .hdrフォーマットで出力 Выходной формат
	char buf[512];
	char s[150];
	char *pointer;
	
	get_time_str(s,end1,startT1);
	
	while((pointer=strchr(s, ':'))!=NULL)
	{
		*pointer='_';
	}
	
	if(mlt_num>1)
	sprintf(buf, "%04d_passes%d_%s-%d.hdr", mutation,mlt_num, s,time(NULL));
	else
	sprintf(buf, "%04d_%s-%d.hdr", mutation, s,time(NULL));
	
	save_hdr_file(buf, image, width, height);
	
	print_time(end1,startT1);
}
