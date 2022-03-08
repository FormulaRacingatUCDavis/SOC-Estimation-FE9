
// This program generates a csv file with data on a simulated state of charge (SOC)
// estimate for a vehicle. It also prints the data in the console (output).
// Time duration of data: 9000s
// Time interval of data: every 10s (this can be modified)
// Created by Ayush Saha, 2/14/2022


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

// time variable declarations
double dt = 1;
int num_seconds = 100; // max duration of sim data
double t = 0;

// circuit model parameters
double Cc = 2400;
double Rc = 0.015;
double Cbat = 18000; // units in Amp * s

// State variables
long double SOCk = 1.0;
long double vk = 4.19;
// next values of state variables
long double SOCkp1;
long double vkp1;

// current load functions
double get_It_constant(double t);
double get_It_toggle(double t);
double get_It_piecewise(double t);
// current load function and value
double (*func_It)(double) = &get_It_piecewise; // choose current function here
double It;

// Linked list struct
struct node {
	double t;
	long double soc;
	long double v;
	double i;
	struct node* next;
};

// WE MUST HAVE A POINTER TO HEAD OF LIST AND POINTER TO CURRENT ELEM
struct node* head = NULL;
struct node* current = NULL;


// INSERT LINK AT THE FIRST LOCATION
void insertFirst(double t_val, long double soc_val, long double v_val, double i_val) {
	//create a link
	struct node* link = (struct node*)malloc(sizeof(struct node));

	link->t = t_val;
	link->soc = soc_val;
	link->v = v_val;
	link->i = i_val;

	current = head = link;
	link->next = NULL;
}

// Insert link at end
struct node* insertAtEnd(double t_val, long double soc_val, long double v_val, double i_val) {
	//create a link
	struct node* link = (struct node*)malloc(sizeof(struct node));
	link->t = t_val;
	link->soc = soc_val;
	link->v = v_val;
	link->i = i_val;
	link->next = NULL;

	current->next = link;
	current = link;
	return link;
}

// PRINTING THE LIST
void printList() {
	struct node* ptr = head;

	while (ptr != NULL) {
		printf("%.4f\t%.9Lf\t%.9Lf\t%.9f\n", ptr->t, ptr->soc, ptr->v, ptr->i);
		ptr = ptr->next;
	}
}

// FINDING LENGTH OF LINKED LIST
int length() {
	int length = 0;
	struct node* current;

	for (current = head; current != NULL; current = current->next) {
		length++;
	}
	return length;
}


// Print the linked list into csv file
int print_list_into_csv(FILE** fp) {
	struct node* ptr = head;
	while (ptr != NULL) {
		fprintf(*fp, "%.9f,%.9Lf,%.9Lf,%.9f\n", ptr->t, ptr->soc, ptr->v, ptr->i);
		ptr = ptr->next;
	}
	return 0;
}

// Discrete solution to
//
//     d/dt(soc) = - It / Cbat
//     d/dt(v) = (It / Cc) - v / (Cc * Rc)
//
// Using finite difference (x_{k+1} - x_{k})/dt approx d/dt(x)
long double func_SOCkp1(long double SOCk, double dt, double It, double Cbat) {
	long double SOCkp1 = SOCk - dt * (It / Cbat);
	return SOCkp1;
}

long double func_vkp1(long double vk, double dt, double It, double Cc, double Rc) {
	long double vkp1 = vk + dt * ((It / Cc) - (vk / (Cc * Rc)));
	return vkp1;
}

// SECTION: options for getting current load: It
// all should have int t as an argument

// Constant
double get_It_constant(double t) {
	return 100;
}

// Toggle on/off
double get_It_toggle(double t) {
	t = (int)t % 10;
	if (t < 5) {
		return 100;
	} else {
		return 0;
	}
}

// Piecewise function
double get_It_piecewise(double t) {
	// First interval
	if ((0 <= t && t < 10) || (30 <= t && t < 40) || (60 <= t && t < 70)) {
		return 300;
	}
	// Second interval
	else if ((10 <= t && t < 20) || (40 <= t && t < 50) || (70 <= t && t < 80)) {
		return 0;
	}
	// Third interval
	else if ((20 <= t && t < 25) || (50 <= t && t < 55) || (80 <= t && t < 85)) {
		return -80;
	}
	// Fourth interval
	else if ((25 <= t && t < 30) || (55 <= t && t < 60) || (85 <= t && t < 100)) {
		return 0;
	}
	return 0;
}

// END SECTION: options for getting current load: It


// DRIVER CODE
int main() {	
	bool start_new_list = false;
	bool done_csv_printing = false;

	// Very first value
	It = (*func_It)(0);
	insertFirst(t, SOCk, vk, It);
	t += dt;

	// Opening the csv file
	FILE* fp;
	fp = fopen("data_sim.csv", "w");
	if (fp == NULL) {
		printf("Can't be opened");
		return 1;
	}

	// THE BIG WHILE LOOP
	while (t <= num_seconds) {

		//----------------GENERATING Y-COOR.----------------

		// get current value at time t
		It = (*func_It)(t);

		// Calculate next SOC and Vc using discrete solution
		SOCkp1 = func_SOCkp1(SOCk, dt, It, Cbat);		
		vkp1 = func_vkp1(vk, dt, It, Cc, Rc);

		// Inserting the (t,soc,v) coordinate into linked list
		if (start_new_list == true) {
			insertFirst(t, SOCkp1, vkp1, It);
			start_new_list = false;
		}
		else {
			insertAtEnd(t, SOCkp1, vkp1, It);
		}
		// Inserting the linked list into the csv
		print_list_into_csv(&fp);
		printList();
		if (t == num_seconds) {
			done_csv_printing = true;
		}
		start_new_list = true;
		current = head;
		// Important variable changes for each run of loop
		vk = vkp1;
		SOCk = SOCkp1;
		t += dt;
	}
	// print last few elements into csv after while loop is over
	if (done_csv_printing == false) {
		print_list_into_csv(&fp);
		printList();
	}

	// Close csv
	fclose(fp);
	fp = NULL;
	return 0;
}
