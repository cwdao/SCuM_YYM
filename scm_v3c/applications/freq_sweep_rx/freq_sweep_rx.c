/**
\brief This program conducts a freqency sweeping test.

After loading this program, SCuM will turn on radio to listen on each
configuration of coarse[5 bits], middle[5 bits] and fine[5 bits] for 15ms.
Meanwhile, an OpenMote will send NUMPKT_PER_CFG packet on one channel every
4ms. SCuM will print out the pkt it received at each config settings.

This program supposes to be run 16 times on for each channel setting on OpenMote
side.
*/

#include <string.h>
#include <math.h>

#include "memory_map.h"
#include "optical.h"
#include "radio.h"
#include "rftimer.h"
#include "scm3c_hw_interface.h"

//=========================== defines =========================================

#define CRC_VALUE (*((unsigned int*)0x0000FFFC))
#define CODE_LENGTH (*((unsigned int*)0x0000FFF8))

#define LENGTH_PACKET 125 + LENGTH_CRC  ///< maximum length is 127 bytes
#define LEN_RX_PKT 20 + LENGTH_CRC      ///< length of rx packet

#define TIMER_PERIOD 500  ///< 500 = 1ms@500kHz

#define NUMPKT_PER_CFG 30
#define STEPS_PER_CONFIG 32

#define NUM_SAMPLES         10

//=========================== variables =======================================

typedef struct {
			uint8_t packet[LENGTH_PACKET];
			uint8_t packet_len;
			int8_t rxpk_rssi;
			uint8_t rxpk_lqi;

			volatile bool rxpk_crc;
			// a flag to mark when to change configure
			volatile bool changeConfig;
			// a flag to avoid change configure during receiving frame
			volatile bool rxFrameStarted;

			volatile uint32_t IF_estimate;
			volatile uint32_t LQI_chip_errors;
			volatile uint32_t cdr_tau_value;

			uint8_t cfg_coarse;
			uint8_t cfg_mid;
			uint8_t cfg_fine;
			
			uint16_t receiveTimes;

            uint32_t 				avg_sample;
            uint32_t 				count_2M;
            uint32_t 				count_LC;
            uint32_t 				count_adc;
            bool 				    freq_setFinish_flag;
            uint8_t 				freq_setTimes_flag;
            bool                    calculate_flag;             //give time for calculating
            uint16_t 				calculate_times_flag;       //give time for calculating
            bool                    RC_count_flag;
            uint8_t			        cb_timer_num;
volatile    uint8_t                 sample_index;
            uint32_t                samples[NUM_SAMPLES];
            bool                    roughEstimateRequireFlag;
            int 				    RC_count;                    //count the number of 2M
            int                     freq_abs;
            //to decrease the time, save the rough estimate parameter
            double                  para_matrix[2][1];
            double                  Inv_equ_c; 
            bool                    statusFlag;    

            uint8_t                 rx_coarse_p1;
            uint8_t                 rx_mid_p1;
            uint8_t                 rx_fine_p1;
            uint8_t                 rx_coarse_p2;
            uint8_t                 rx_mid_p2;
            uint8_t                 rx_fine_p2;
} app_vars_t;

app_vars_t app_vars;

double    freqRXTargetList[16] = {1.251 ,1.254 ,1.257 ,1.259 ,1.262 ,1.264 ,1.267 ,1.270 ,1.272 ,1.275 ,1.277 ,1.280 ,1.283 ,1.285 ,1.288 ,1.290};

//=========================== prototypes ======================================

void        cb_startFrame_rx(uint32_t timestamp);
void        cb_endFrame_rx(uint32_t timestamp);
void        cb_timer(void);

void        __PresiseEstimate(uint16_t channelTarget);
void        course_estimate(uint16_t channelTarget);
void        Gaussian_elimination(int x_matrix_raw[2][2], double x_matrix_inverse[2][2]);
void        matrixMultiplication(double x_matrix_inverse[][2], int y_matrix[][1], double para_matrix[][1]);
void        Solve_Equation(double para_matrix[][1], double Inv_equ_c, double x_matrix_inverse[][2]);
uint32_t    average_sample(void);
//=========================== main ============================================

int main(void) {
    uint32_t calc_crc;

    uint8_t i;
    uint8_t j;
    uint8_t offset;
    uint16_t channelTarget;

    memset(&app_vars, 0, sizeof(app_vars_t));

    printf("Initializing...");

    // Set up mote configuration
    // This function handles all the analog scan chain setup
    initialize_mote();

    // radio_setStartFrameRxCb(cb_startFrame_rx);
    // radio_setEndFrameRxCb(cb_endFrame_rx);
    rftimer_set_callback(cb_timer);

    // Disable interrupts for the radio and rftimer
    radio_disable_interrupts();
    rftimer_disable_interrupts();

    // Check CRC to ensure there were no errors during optical programming
    printf("\r\n-------------------\r\n");
    printf("Validating program integrity...");

    calc_crc = crc32c(0x0000, CODE_LENGTH);

    if (calc_crc == CRC_VALUE) {
        printf("CRC OK\r\n");
    } else {
        printf(
            "\r\nProgramming Error - CRC DOES NOT MATCH - Halting "
            "Execution\r\n");
        while (1)
            ;
    }

    // Debug output
    // printf("\r\nCode length is %u bytes",code_length);
    // printf("\r\nCRC calculated by SCM is: 0x%X",calc_crc);

    // printf("done\r\n");

    // After bootloading the next thing that happens is frequency calibration
    // using optical
    printf("Calibrating frequencies...\r\n");

    // Initial frequency calibration will tune the frequencies for HCLK, the
    // RX/TX chip clocks, and the LO


    // For the LO, calibration for RX channel 11, so turn on AUX, IF, and LO
    // LDOs by calling radio rxEnable
    radio_rxEnable();

    // Enable optical SFD interrupt for optical calibration
    optical_enable();

    // Wait for optical cal to finish
    optical_enableLCCalibration();
    while (optical_getCalibrationFinished() == 0) ;
    //default as tx mode
    // Turn on LO, DIV, PA, and IF
    ANALOG_CFG_REG__10 = 0x78;

    printf("Cal complete\r\n");

    // Enable interrupts for the radio FSM
    radio_enable_interrupts();
    channelTarget = 11;
    // configure
    course_estimate(channelTarget);
    __PresiseEstimate(channelTarget);

    while (1) {
        // loop through all configuration
        // for (app_vars.cfg_coarse = 24; app_vars.cfg_coarse < 27;
        //      app_vars.cfg_coarse++) {
        //     for (app_vars.cfg_mid = 15; app_vars.cfg_mid < STEPS_PER_CONFIG;
        //          app_vars.cfg_mid++) {
        //         for (app_vars.cfg_fine = 0;
        //              app_vars.cfg_fine < STEPS_PER_CONFIG;
        //              app_vars.cfg_fine++) {
        //             printf("coarse=%d, middle=%d, fine=%d\r\n", app_vars.cfg_coarse,app_vars.cfg_mid,app_vars.cfg_fine);
		// 							  app_vars.receiveTimes = 0;
        //             for (i = 0; i < NUMPKT_PER_CFG; i++) {
        //                 while (app_vars.rxFrameStarted == true)
        //                     ;
        //                 radio_rfOff();
        //                 LC_FREQCHANGE(app_vars.cfg_coarse, app_vars.cfg_mid,
        //                               app_vars.cfg_fine);
        //                 radio_rxEnable();
        //                 radio_rxNow();
        //                 rftimer_setCompareIn(rftimer_readCounter() +
        //                                      TIMER_PERIOD);
        //                 app_vars.changeConfig = false;
        //                 while (app_vars.changeConfig == false)
        //                     ;
												
        //             }
		// 								if(app_vars.receiveTimes != 0){
		// 								printf("received%d\r\n",app_vars.receiveTimes);
		// 								}
        //         }
        //     }
        // }
        
    }
}

//=========================== public ==========================================

//=========================== private =========================================

void cb_startFrame_rx(uint32_t timestamp) { app_vars.rxFrameStarted = true;}

void cb_endFrame_rx(uint32_t timestamp) {
    uint8_t i;

    radio_getReceivedFrame(&(app_vars.packet[0]), &app_vars.packet_len,
                           sizeof(app_vars.packet), &app_vars.rxpk_rssi,
                           &app_vars.rxpk_lqi);

    radio_rfOff();

    if (app_vars.packet_len == LEN_RX_PKT && (radio_getCrcOk())) {
        // Only record IF estimate, LQI, and CDR tau for valid packets
        app_vars.IF_estimate = radio_getIFestimate();
        app_vars.LQI_chip_errors = radio_getLQIchipErrors();

//        printf("pkt received on ch%d %c%c%c%c.%d.%d.%d\r\n", app_vars.packet[4],
//               app_vars.packet[0], app_vars.packet[1], app_vars.packet[2],
//               app_vars.packet[3], app_vars.cfg_coarse, app_vars.cfg_mid,
//               app_vars.cfg_fine);
				app_vars.receiveTimes++;
        app_vars.packet_len = 0;
        memset(&app_vars.packet[0], 0, LENGTH_PACKET);
    }

    radio_rxEnable();
    radio_rxNow();

    app_vars.rxFrameStarted = false;
}

void cb_timer(void) {   
    if(app_vars.statusFlag == 0){
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        read_counters_3B(&app_vars.count_2M,&app_vars.count_LC,&app_vars.count_adc);
        if (app_vars.RC_count_flag==0)
        {
            app_vars.samples[app_vars.sample_index] = app_vars.count_2M;
            app_vars.sample_index++;
            if (app_vars.sample_index==NUM_SAMPLES) {
                app_vars.sample_index = 0;
                app_vars.RC_count = average_sample();
                app_vars.cb_timer_num++;
            }
            if (app_vars.cb_timer_num==4)
            {
                app_vars.RC_count_flag = 1;
            }
            
            // printf("here");
        }
        if (app_vars.RC_count_flag==1){
            app_vars.samples[app_vars.sample_index] = app_vars.count_LC;
            app_vars.sample_index++;
            if (app_vars.sample_index==NUM_SAMPLES) {
                app_vars.sample_index = 0;
                app_vars.avg_sample = average_sample();
        //			printf("%d \r\n",app_vars.avg_sample);
                    
                if(app_vars.freq_setFinish_flag == 0 && app_vars.freq_setTimes_flag == 5){
                    app_vars.freq_setFinish_flag = 1;
                    app_vars.freq_setTimes_flag = 0;
                }
                app_vars.freq_setTimes_flag++;
            }
            if (app_vars.calculate_flag == 1)
            {
                app_vars.calculate_times_flag++;
                if (app_vars.calculate_times_flag==30)
                {
                    app_vars.calculate_times_flag=0;
                    app_vars.calculate_flag = 0;
                }
                
            }
        }
    }
    else if (app_vars.statusFlag == 1)
    {
        app_vars.changeConfig = true; 
    }
}


void course_estimate(uint16_t channelTarget){
   
    //用无符号数总会出现奇怪的问题
    int p1L_count, p1R_count, p2L_count, p2R_count, p3L_count, p3R_count;
    int p1x, p1y, p2x, p2y, p3x, p3y;
    int x_matrix[2][2];
    int y_matrix[2][1];
    double x_matrix_inverse[2][2];
    double para_matrix[2][1];
    double Inv_equ_c;

    if(app_vars.roughEstimateRequireFlag == 1){
        p1x, p1y, p2x, p2y, p3x, p3y = 0;
        Inv_equ_c = 0;
    }
    // double test1, test2, test3;
    rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
    while (app_vars.RC_count_flag==0);
    
    LC_FREQCHANGE(1,31,31);
    rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
    while(app_vars.freq_setFinish_flag == 0);
    p1L_count = app_vars.avg_sample;
    while (p1L_count <= 1000 || p1L_count >= 3000)       
    {
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        while(app_vars.freq_setFinish_flag == 0);
        printf("cal error, waiting..");
        p1L_count = app_vars.avg_sample;
    }
    
    printf("p1L_count=%d\r\n",p1L_count);

    LC_FREQCHANGE(2,0,0);
    rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
    app_vars.freq_setFinish_flag = 0;
    while(app_vars.freq_setFinish_flag == 0);
    p1R_count = app_vars.avg_sample;
    // printf("p1R_count=%d\r\n",p1R_count);
    
    LC_FREQCHANGE(15,31,31);
    rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
    app_vars.freq_setFinish_flag = 0;
    while(app_vars.freq_setFinish_flag == 0);
    p2L_count = app_vars.avg_sample;
    // printf("p2L_count=%d\r\n",p2L_count);
    
    LC_FREQCHANGE(16,0,0);
    rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
    app_vars.freq_setFinish_flag = 0;
    while(app_vars.freq_setFinish_flag == 0);
    p2R_count = app_vars.avg_sample;
    // printf("p2R_count=%d\r\n",p2R_count);
    
    LC_FREQCHANGE(30,31,31);
    rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
    app_vars.freq_setFinish_flag = 0;
    while(app_vars.freq_setFinish_flag == 0);
    p3L_count = app_vars.avg_sample;
    // printf("p3L_count=%d\r\n",p3L_count);
    
    LC_FREQCHANGE(31,0,0);
    rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
    app_vars.freq_setFinish_flag = 0;
    while(app_vars.freq_setFinish_flag == 0);
    p3R_count = app_vars.avg_sample;
    // printf("RC_count=%d\r\n",p3R_count);
    // Inv_equ_c = (double)app_vars.RC_count;
    // // printf("ah=%f",Inv_equ_c);
    printf("RC=%d\r\n",app_vars.RC_count);

    p1y = (p1L_count + p1R_count)/2;
    p2y = (p2L_count + p2R_count)/2;
    p3y = (p3L_count + p3R_count)/2;
    p1x = 2;
    p2x = 16;
    p3x = 31;

    // 修改 x_matrix
    x_matrix[0][0] = ((p2x * p2x) - (p1x * p1x));
    x_matrix[0][1] = (p2x - p1x);
    x_matrix[1][0] = ((p3x * p3x) - (p1x * p1x));
    x_matrix[1][1] = (p3x - p1x);
    // printf("xmatrix = %d,%d,%d,%d\r\n", x_matrix[0][0], x_matrix[0][1], x_matrix[1][0], x_matrix[1][1]);
    //0427 浮点运算算力太差，且使用标号为x会导致结果太大
    //能用整数尽可能用整数
    //尝试直接以x=2, x=16, x=31进行计算
    //反解出来只需要和.5比较就可以

    //0428 部分板子存在course code小的时候LC不准确的情况
    //尝试使用更大的course code来实现拟合

    // 修改 y_matrix
    y_matrix[0][0] = p2y - p1y;
    y_matrix[1][0] = p3y - p1y;
    // printf("y_matrix = %d,%d\r\n", y_matrix[0][0], y_matrix[1][0]);
    // printf("xmatrix = %d,%d,%d,%d\r\n", x_matrix[0][0], x_matrix[0][1], x_matrix[1][0], x_matrix[1][1]);
    Gaussian_elimination(x_matrix, x_matrix_inverse);
    matrixMultiplication(x_matrix_inverse, y_matrix, para_matrix);
    // printf("para_matrix = %f,%f\r\n", para_matrix[0][0], para_matrix[1][0]);
    //x_matrix_inverse不再使用，帮助计算存储中间值
    x_matrix_inverse[0][0] = para_matrix[0][0] * (p2x *p2x);
    x_matrix_inverse[0][1] = para_matrix[1][0] * p2x;
    x_matrix_inverse[1][0] = x_matrix_inverse[0][0] + x_matrix_inverse[0][1];
    Inv_equ_c = p2y - x_matrix_inverse[1][0];
    printf("Inv_equ_c=%f\r\n",Inv_equ_c);
    app_vars.Inv_equ_c = Inv_equ_c;
    app_vars.para_matrix[0][0] = para_matrix[0][0];
    app_vars.para_matrix[1][0] = para_matrix[1][0];
    y_matrix[0][0] = 0;
    y_matrix[1][0] = 0;
}

void __PresiseEstimate(uint16_t channelTarget){
    int p1L_count, p1R_count, p2L_count, p2R_count, p3L_count, p3R_count;
    int p1x, p1y, p2x, p2y, p3x, p3y;
    double x_matrix_inverse[2][2];
    double assistMatrix[2][1];  //assist calculating, sometime directive calculating will stuck
    
    printf("begin precise stimate\r\n");
    app_vars.freq_abs = (int)(freqRXTargetList[channelTarget-11]*app_vars.RC_count);
    printf("freq=%d\r\n",app_vars.freq_abs);
    Solve_Equation(app_vars.para_matrix, app_vars.Inv_equ_c,x_matrix_inverse);
    //开始分区，p1L_count, p1R_count, p2L_count, p2R_count, p3L_count, p3R_count释放，帮助计算
    if ((int)x_matrix_inverse[0][0]%10 <=5)
    {
        p1R_count = (int)x_matrix_inverse[1][0];
        printf("p1r=%d\r\n",p1R_count);
        p1L_count = p1R_count - 1;
        LC_FREQCHANGE(p1L_count,16,0);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p1x = app_vars.avg_sample;

        LC_FREQCHANGE(p1L_count,31,31);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p1y = app_vars.avg_sample;

        LC_FREQCHANGE(p1R_count,0,0);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p2x = app_vars.avg_sample;

        LC_FREQCHANGE(p1R_count,15,31);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p2y = app_vars.avg_sample;
        app_vars.rx_coarse_p1 = p1L_count;
        app_vars.rx_coarse_p2 = p1R_count;
    }
    else if ((int)x_matrix_inverse[0][0]%10 >5)
    {
        p1L_count = (int)x_matrix_inverse[1][0];
        printf("p1r=%d\r\n",p1R_count);
        p1R_count = p1L_count + 1;
        LC_FREQCHANGE(p1L_count,16,0);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p1x = app_vars.avg_sample;

        LC_FREQCHANGE(p1L_count,31,31);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p1y = app_vars.avg_sample;

        LC_FREQCHANGE(p1R_count,0,0);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p2x = app_vars.avg_sample;

        LC_FREQCHANGE(p1R_count,15,31);
        rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
        app_vars.freq_setFinish_flag = 0;
        while(app_vars.freq_setFinish_flag == 0);
        p2y = app_vars.avg_sample;
        app_vars.rx_coarse_p1 = p1L_count;
        app_vars.rx_coarse_p2 = p1R_count;
    }
    //存在隐藏问题，如果比两个都小，待修正
    //可以优化，如果比其中一个大，顺延到下一个course的下区
    while ((app_vars.freq_abs > p1y && app_vars.freq_abs > p2y))
    {
        x_matrix_inverse[1][0] = x_matrix_inverse[1][0] + 1;
        if ((int)x_matrix_inverse[0][0]%10 <=5)
        {
            p1R_count = (int)x_matrix_inverse[1][0];
            printf("p1r=%d\r\n",p1R_count);
            p1L_count = p1R_count - 1;
            LC_FREQCHANGE(p1L_count,16,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1x = app_vars.avg_sample;

            LC_FREQCHANGE(p1L_count,31,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1y = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,0,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2x = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,15,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2y = app_vars.avg_sample;
            app_vars.rx_coarse_p1 = p1L_count;
            app_vars.rx_coarse_p2 = p1R_count;
        }
        else if ((int)x_matrix_inverse[0][0]%10 >5)
        {
            p1L_count = (int)x_matrix_inverse[1][0];
            printf("p1r=%d\r\n",p1R_count);
            p1R_count = p1L_count + 1;
            LC_FREQCHANGE(p1L_count,16,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1x = app_vars.avg_sample;

            LC_FREQCHANGE(p1L_count,31,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1y = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,0,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2x = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,15,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2y = app_vars.avg_sample;

            app_vars.rx_coarse_p1 = p1L_count;
            app_vars.rx_coarse_p2 = p1R_count;
        }
        printf("p1x=%d, p1y=%d, p2x=%d, p2y=%d\r\n",p1x, p1y, p2x, p2y);
    }
    //如果比任意一个小
    while ((app_vars.freq_abs < p1x && app_vars.freq_abs < p2x))
    {
        x_matrix_inverse[1][0] = x_matrix_inverse[1][0] - 1;
        if ((int)x_matrix_inverse[0][0]%10 <=5)
        {
            p1R_count = (int)x_matrix_inverse[1][0];
            printf("p1r=%d\r\n",p1R_count);
            p1L_count = p1R_count - 1;
            LC_FREQCHANGE(p1L_count,16,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1x = app_vars.avg_sample;

            LC_FREQCHANGE(p1L_count,31,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1y = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,0,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2x = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,15,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2y = app_vars.avg_sample;
            app_vars.rx_coarse_p1 = p1L_count;
            app_vars.rx_coarse_p2 = p1R_count;
        }
        else if ((int)x_matrix_inverse[0][0]%10 >5)
        {
            p1L_count = (int)x_matrix_inverse[1][0];
            printf("p1r=%d\r\n",p1R_count);
            p1R_count = p1L_count + 1;
            LC_FREQCHANGE(p1L_count,16,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1x = app_vars.avg_sample;

            LC_FREQCHANGE(p1L_count,31,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p1y = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,0,0);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2x = app_vars.avg_sample;

            LC_FREQCHANGE(p1R_count,15,31);
            rftimer_setCompareIn(rftimer_readCounter()+TIMER_PERIOD);
            app_vars.freq_setFinish_flag = 0;
            while(app_vars.freq_setFinish_flag == 0);
            p2y = app_vars.avg_sample;

            app_vars.rx_coarse_p1 = p1L_count;
            app_vars.rx_coarse_p2 = p1R_count;
        }
        printf("p1x=%d, p1y=%d, p2x=%d, p2y=%d\r\n",p1x, p1y, p2x, p2y);
    }
    
    //打印获取的数据
    printf("p1x=%d, p1y=%d, p2x=%d, p2y=%d\r\n",p1x, p1y, p2x, p2y);
    //都围绕在course的半区内，p1区使用16-32，p2区使用0-16
    // p1区
    p2L_count = p1y - p1x;
    p2R_count = 32 - 16;
    x_matrix_inverse[0][0] = (double)p2L_count/p2R_count; //k
    x_matrix_inverse[0][1] = x_matrix_inverse[0][0]*16; //k*x
    x_matrix_inverse[1][0] = p1x - x_matrix_inverse[0][1]; //b=y-k*x
    x_matrix_inverse[0][1] = app_vars.freq_abs - x_matrix_inverse[1][0]; //y-b
    assistMatrix[0][0] = x_matrix_inverse[0][1]/x_matrix_inverse[0][0]; //p1 set = (y-b)/k

    // p2区
    p3L_count = p2y - p2x;
    p3R_count = 16 - 0;
    x_matrix_inverse[0][0] = (double)p3L_count/p3R_count; //k
    x_matrix_inverse[0][1] = x_matrix_inverse[0][0]*16; //k*x
    x_matrix_inverse[1][0] = p2y - x_matrix_inverse[0][1];// b = y-kx
    x_matrix_inverse[0][1] = app_vars.freq_abs - x_matrix_inverse[1][0]; //y-b
    assistMatrix[1][0] = x_matrix_inverse[0][1]/x_matrix_inverse[0][0]; //p2 set = (y-b)/k

    //下面需要验证是否在该范围中
    //可能遇到的问题是LC_count与发射频率的误差

    //输出结果
    if ((int)assistMatrix[0][0]>=16 && (int)assistMatrix[0][0]<=32)
    {
        p1x = (int)assistMatrix[0][0];
        x_matrix_inverse[0][1] = assistMatrix[0][0] - p1x;
        p1y = (int)(x_matrix_inverse[0][1]*31);
    }
    else
    {
        p1x = 0;
        p1y = 0;
    }
    if ((int)assistMatrix[1][0]>=0 && (int)assistMatrix[1][0]<=16)
    {
        p2x = (int)assistMatrix[1][0];
        x_matrix_inverse[1][1] = assistMatrix[1][0] - p2x;
        p2y = (int)(x_matrix_inverse[1][1]*31);
    }
    else
    {
        p2x = 0;
        p2y = 0;
    }
    printf("set1=%d.%d.%d\r\n",app_vars.rx_coarse_p1,p1x,p1y);
    printf("set2=%d.%d.%d\r\n",app_vars.rx_coarse_p2,p2x,p2y); 
    app_vars.rx_mid_p1 = p1x;
    app_vars.rx_fine_p1 = p1y;
    app_vars.rx_mid_p2 = p2x;
    app_vars.rx_fine_p2 = p2y;
    app_vars.statusFlag = 1;
    
}


void Gaussian_elimination(int x_matrix_raw[2][2], double x_matrix_inverse[2][2]) {
    // 计算矩阵的行列式
    double determinant;

    determinant = x_matrix_raw[0][0] * x_matrix_raw[1][1] - x_matrix_raw[0][1] * x_matrix_raw[1][0];
    
    // 检查行列式是否为0，若为0则矩阵不可逆
    if (determinant == 0) {
        printf("Matrix is not invertible.\n");
        return;
    }

    // 计算逆矩阵的每个元素
    x_matrix_inverse[0][0] = x_matrix_raw[1][1] / determinant;
    x_matrix_inverse[0][1] = -x_matrix_raw[0][1] / determinant;
    x_matrix_inverse[1][0] = -x_matrix_raw[1][0] / determinant;
    x_matrix_inverse[1][1] = x_matrix_raw[0][0] / determinant;
    app_vars.calculate_flag = 1;
    while(app_vars.calculate_flag == 1){
    }; 
}

/**************************/
//Function: multiple two matrix
//Parameter: x_matrix_raw[][],  x_matrix_inversep[][]
//Return: None (return matrix directively)
//Call: course_estimate 
//Date: 31.May.2024
/*************************/
void matrixMultiplication(double x_matrix_inverse[][2], int y_matrix[][1], double para_matrix[][1]) {
    int i, j, k;
    
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 1; j++) {
            para_matrix[i][j] = 0;
            for (k = 0; k < 2; k++) {
                para_matrix[i][j] += x_matrix_inverse[i][k] * y_matrix[k][j];
            }
        }
    }
}

void Solve_Equation(double para_matrix[][1], double Inv_equ_c, double x_matrix_inverse[][2]) {
    Inv_equ_c = Inv_equ_c - app_vars.freq_abs;
    //同样使用x_matrix_inverse帮助计算存储中间值，以减少内存
    x_matrix_inverse[0][0] = para_matrix[1][0]*para_matrix[1][0];
    x_matrix_inverse[0][1] = 4 * para_matrix[0][0]*Inv_equ_c;
    x_matrix_inverse[1][0] = x_matrix_inverse[0][0] - x_matrix_inverse[0][1];
    x_matrix_inverse[1][1] = sqrt(x_matrix_inverse[1][0]);

    //x_matrix_inverse[0][0],[0][1],[1][0]重新释放 
    x_matrix_inverse[0][0] = 2*para_matrix[0][0];
    x_matrix_inverse[0][1] = 1/x_matrix_inverse[0][0];
    x_matrix_inverse[1][0] = x_matrix_inverse[0][1] * (-para_matrix[1][0]+x_matrix_inverse[1][1]);
    //只需要x大于0的根
    //在此处直接乘10，看除10的余数大于5还是小于5
    x_matrix_inverse[0][0] = x_matrix_inverse[1][0] * 10;
    // printf("solve = %f\r\n", x_matrix_inverse[1][0]);
}

uint32_t  average_sample(void){
    uint8_t i;
    uint32_t avg;
    
    avg = 0;
    for (i=0;i<NUM_SAMPLES;i++) {
        avg += app_vars.samples[i];
    }
    avg = avg/NUM_SAMPLES;
    return avg;
}