#include "pch.h"

#include <stdio.h>
#include <Windows.h>
#include <GL/GL.h>
#include <GL/glut.h>

#include <HD/hd.h>
#include <HDU/hduError.h>
#include <HDU/hduMatrix.h>
#include <math.h>

const double mass = 1.0;
const double gravity = 9.80665;
const double dt = 0.1;
const double k = 10.0;
const double natural_length = 0.0;
const double first_accel = gravity;
const double first_velocity = 0.0;
const double first_position = 0.0;

double R_x = first_position;
double R_v = first_velocity;
double t = 0;
double accel = first_accel;
double b = 0.1;

GLfloat object_color[4] = { 0.0, 1.0, 0.0, 1.0 };
GLfloat line_color[4] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_color[4] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_position[4] = { 100.0, 100.0, 100.0, 1.0 };

// --- グローバル変数 ---
hduVector3Dd gCenterOfStylusSphere, gCenterOfGodSphere, gForce;
HHD ghHD = HD_INVALID_HANDLE;
HDSchedulerHandle gSchedulerCallback = HD_INVALID_HANDLE;

// --- 触覚の空間をグラフィックスの空間と同じにする ---
void resize(int w, int h)
{
	glMatrixMode(GL_PROJECTION); // 投影投影モードへ
	glLoadIdentity(); // 投影変換の変換行列を単位行列で初期化
	glOrtho(-100.0, 100.0, -100.0, 100.0, -100.0, 100.0); //各軸指定した範囲で囲まれる立方体の範囲を平行投影
}

// --- グラフィックスパイプライン，光源の初期化 ---
void doGraphicsState()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // カラーバッファとデプスバッファをクリア
	glEnable(GL_COLOR_MATERIAL); // 物体色の有効化
	glEnable(GL_LIGHTING); // ライティング（照光）処理を有効化
	glEnable(GL_NORMALIZE); // 法線ベクトルを正規化（単位ベクトル化）
	glShadeModel(GL_SMOOTH); // スムースシェーディングを有効化
	GLfloat lightZeroPosition[] = { 10.0f, 4.0f, 100.0f, 0.0f }; // 光源０の位置
	GLfloat lightZeroColor[] = { 0.6f, 0.6f, 0.6f, 1.0f }; // 光源０の色
	GLfloat lightOnePosition[] = { -1.0f, -2.0f, -100.0f, 0.0f }; // 光源１の位置
	GLfloat lightOneColor[] = { 0.6f, 0.6f, 0.6f, 1.0f }; // 光源１の色
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE); // より正確な陰影付けを行う
	glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
	//glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
	//glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
	glEnable(GL_LIGHT0); // 光源０の有効化
	//glEnable(GL_LIGHT1); // 光源１の有効化
	glEnable(GL_DEPTH_TEST);
}

// --- GLUTライブラリから定期的に呼ばれる関数 ---
void idle(void)
{
	glutPostRedisplay(); // グラフィックスの再描画
	if (!hdWaitForCompletion(gSchedulerCallback, HD_WAIT_CHECK_STATUS)) {
		printf("メインスケジューラが終了しました。\n何かキーを押すと終了します。\n");
		getchar();
		exit(-1);
	}
}

double F1(double t, double x, double v) {
	return -k / mass * x - b * v / mass;
}

double F2(double t, double x, double v)
{
	return v;
}

void Runge_Kutta() {
	t += dt;

	double k1_x = F2(t, R_x, R_v);
	double k1_v = F1(t, R_x, R_v);
	double k2_x = F2(t + dt / 2, R_x + k1_x * dt / 2, R_v + k1_v * dt / 2);
	double k2_v = F1(t + dt / 2, R_x + k1_x * dt / 2, R_v + k1_v * dt / 2);
	double k3_x = F2(t + dt / 2, R_x + k2_x * dt / 2, R_v + k2_v * dt / 2);
	double k3_v = F1(t + dt / 2, R_x + k2_x * dt / 2, R_v + k2_v * dt / 2);
	double k4_x = F2(t + dt, R_x + k3_x * dt, R_v + k3_v * dt / 2);
	double k4_v = F1(t + dt, R_x + k3_x * dt, R_v + k3_v * dt / 2);

	R_v += dt / 6 * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
	R_x += dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);

	printf("psition: %f\t velocity: %f\n", R_x, R_v);
}

#define SSR 10.0
#define FSR 60.0
#define OBJECT_MASS 30.0
#define STIFFNESS 0.10
#define EPSILON 0.00001

// --- GodObjectの座標を更新 ---
void updateEffectorPosition(void)
{
	double m_currentDistance;
	hduVector3Dd m_centerToEffector = gCenterOfStylusSphere;
	m_currentDistance = sqrt(m_centerToEffector[0] * m_centerToEffector[0] + ((m_centerToEffector[1]+R_x) * (m_centerToEffector[1]+R_x)) + m_centerToEffector[2] * m_centerToEffector[2]);

	/*if (m_currentDistance > SSR + FSR) {
		gCenterOfGodSphere = gCenterOfStylusSphere;
		gForce.set(0.0, 0.0, 0.0);
	}
	else {
		if (m_currentDistance > EPSILON) {
			double scale = (SSR + FSR) / m_currentDistance;
			gCenterOfGodSphere = scale * m_centerToEffector;
			gForce = STIFFNESS * (gCenterOfGodSphere - gCenterOfStylusSphere);

		}
	}*/

	if (m_currentDistance > SSR + SSR) {
		gCenterOfGodSphere = gCenterOfStylusSphere;
		gForce.set(0.0, 0.0, 0.0);
	}
	else {
		if (m_currentDistance > EPSILON) {
			double scale = (SSR + SSR) / m_currentDistance;
			gCenterOfGodSphere = scale * m_centerToEffector;
			gForce = STIFFNESS * (gCenterOfGodSphere - gCenterOfStylusSphere);
			R_x += gForce[1];
		}
	}
}

// --- 力覚スケジューラ ---
HDCallbackCode HDCALLBACK ContactCB(void* data)
{


	HHD hHD = hdGetCurrentDevice();
	hdBeginFrame(hHD);
	hdGetDoublev(HD_CURRENT_POSITION, gCenterOfStylusSphere);
	updateEffectorPosition();
	hdSetDoublev(HD_CURRENT_FORCE, gForce);
	hdEndFrame(hHD);
	HDErrorInfo error;
	if (HD_DEVICE_ERROR(error = hdGetError())) {
		hduPrintError(stderr, &error, "力覚スケジューラ内でエラーが発生しました。");
		if (hduIsSchedulerError(&error)) {
			return HD_CALLBACK_DONE;
		}
	}
	return HD_CALLBACK_CONTINUE;
}

// グラフィックス（球体）の表示
void display()
{
	doGraphicsState();
	glDisable(GL_CULL_FACE); // カリングの無効化

	glPushMatrix();
	
	Runge_Kutta();

	glPushMatrix();
	glColor4fv(line_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, line_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, line_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, line_color);
	glTranslated(0.0, 50.0, 0.0);
	glutSolidSphere(SSR/2, 20, 20);
	glPopMatrix();

	glLineWidth(10.0);
	glBegin(GL_LINES);
	glColor4fv(line_color);
	glVertex3d(0.0, 50.0, 0.0);
	glVertex3d(0.0, R_x, 0.0);
	glEnd();

	glPushMatrix();
	glColor4fv(object_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, object_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, object_color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, object_color);
	glTranslated(0.0, R_x, 0.0);
	glutSolidSphere(SSR, 20, 20);
	glPopMatrix();;

	glPushMatrix();
	glTranslated(gCenterOfGodSphere[0], gCenterOfGodSphere[1], gCenterOfGodSphere[2]);
	glColor4d(1, 1, 1, 1.0);
	glutSolidSphere(SSR, 20, 20); // 物体表面にとどまる球(GodObject)
	glPopMatrix();

	//The device-controlled sphere.
	glPushMatrix();
	glTranslated(gCenterOfStylusSphere[0], gCenterOfStylusSphere[1], gCenterOfStylusSphere[2]);
	glColor4d(1.0, 1.0, 0.0, 1.0);
	glutWireSphere(SSR, 20, 20); // 物体内部に侵入する球
	glPopMatrix();


	glutSwapBuffers();
}
// --- 終了処理 ---
void exitHandler()
{
	hdStopScheduler();
	if (ghHD != HD_INVALID_HANDLE) {
		hdDisableDevice(ghHD);
		ghHD = HD_INVALID_HANDLE;
	}
	printf("終了。\n");
}
// --- キーボード入力時のコールバック関数 ----
void keyboard(unsigned char key, int x, int y)
{
	if (key == 'q') exit(0); // 'q'キーが押されたらプログラムを終了
}
// --- メイン関数 ---
int main(int argc, char* argv[])
{

	HDErrorInfo error;
	printf("起動します\n");
	atexit(exitHandler);
	ghHD = hdInitDevice(HD_DEFAULT_DEVICE);
	hdEnable(HD_FORCE_OUTPUT);
	hdEnable(HD_MAX_FORCE_CLAMPING);
	hdStartScheduler();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
	glutCreateWindow("hellohaptics");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	gSchedulerCallback = hdScheduleAsynchronous(ContactCB, NULL, HD_DEFAULT_SCHEDULER_PRIORITY);
	glutMainLoop(); // GLUTのメインループ
	return 0;
}